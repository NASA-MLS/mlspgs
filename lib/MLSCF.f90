! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSCF                    ! MLSCF stuff including reading
!=============================================================================

USE SDPToolkit
USE Hdf
USE MLSMessageModule
USE MLSStrings
USE MLSCommon

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
INTEGER, PARAMETER :: MaxNOMlscfKeys = 38
INTEGER, PARAMETER :: MaxCharValueLen =80
INTEGER, PARAMETER :: UnitsLen = 10
INTEGER, PARAMETER :: MlscfEntryLen = 25
INTEGER, PARAMETER :: MaxNoKeysPerEntry = 20
INTEGER, PARAMETER :: MaxNoEntriesPerSection = 26

TYPE mlscfKey_T
  CHARACTER (len=MaxKeyLen) :: Keyword
  CHARACTER (len=MaxTypeLen):: TYPE
  INTEGER :: Keylen
END TYPE mlscfKey_T



TYPE MlscfCell_T

   ! This datatype defines a Mlscf "cell", i.e. a contruction of type:
   ! "keyword=value" or "keyword=LowerBound..UpperBound[units]"

   ! First the "keyword" holder:

   CHARACTER (len=MaxKeyLen) :: Keyword
   ! Value holder for "real" type

   REAL (R8) :: RealValue

   ! Value holder for integer type

   INTEGER :: IntValue

   ! Value holder for "character" type

   CHARACTER (len= MaxCharValueLen) :: CharValue

   ! Value holder for "units" field

   CHARACTER (len=UnitsLen) :: Units

   ! Value holder for upper and lower bounds fields

   REAL (R8) :: RangeLowerBound
   REAL (R8) :: RangeUpperBound


END TYPE MlscfCell_T


TYPE MlscfEntry_T

   ! This datatype defines am Mlscf entry, i.e. a name followed by
   ! a list of "keyword=value" cells
 
   ! Entry name, e.g. aprioriTemp

   CHARACTER (len=MlscfEntryLen) MlscfEntryName

   ! Number of "keyword=value" cells

   INTEGER MlscfEntryNoKeys

   !Cells holder 

   TYPE (MlscfCell_T), DIMENSION(MaxNoKeysPerEntry) :: Cells

END TYPE MlscfEntry_T


TYPE MlscfSection_T


   ! For each section:
   ! Name, Number of entries within section
   CHARACTER (len=MlscfEntryLen) MlscfSectionName
   INTEGER NoSectionEntries

   ! The actual entries:

   TYPE (MlscfEntry_T), DIMENSION(MaxNoEntriesPerSection) :: Entries ! (NoSectionEntries)


END TYPE MlscfSection_T




TYPE Mlscf_T

! This datatype defines the mlscf as an array of MlscfSections

  INTEGER NoSections

  TYPE (MlscfSection_T), POINTER, DIMENSION (:) :: Sections ! (NoSections) 

END TYPE mlscf_T



! --------------------------------------------------------------------------

CONTAINS


! This function extracts a named section from the mlscf

!========================================
FUNCTION GetMLSCFSection(mlscfInfo,sectionName)
!========================================

  ! Dummy arguments
  TYPE (mlscf_T), INTENT(IN) :: mlscfInfo
  CHARACTER (LEN=*), INTENT(IN) :: sectionName

  ! Function result
  TYPE (mlscfSection_T) :: GetMLSCFSection

  ! Local variables
  INTEGER :: sindex
  LOGICAL :: found

  ! Executable code
  sindex = 1
  found = .false.

  do while (.not. found .and. sindex <= mlscfInfo%NoSections)
     if(sectionName(1:) == mlscfInfo%sections(sindex)%mlscfSectionName(1:))then
        found = .true.
     else
        sindex = sindex + 1
     end if
  end do 

  IF (.not. found) CALL MLSMessage(MLSMSG_Error,ModuleName,&
       & "No such mlscf section:"//sectionName)

  GetMLSCFSection=mlscfInfo%sections(sindex)
END FUNCTION GetMLSCFSection
  
!========================================
SUBROUTINE Strip_Blanks (line, nc, ncnb)
!========================================
 

CHARACTER (LEN=*) :: line
INTEGER :: nc, ncnb
INTEGER :: i, ncc

LOGICAL :: done = .FALSE.

done = .FALSE.
ncc = nc

! check for blank line
line = TRIM(line(1:ncc))
ncnb = LEN_TRIM(line)
IF (ncnb ==  0)RETURN

! ignore trailing comments

i= INDEX(line, ';')
IF(  i /= 0)THEN
  ncc = i-1
  ncnb = ncc

! check again for blank line
  line = TRIM(line(1:ncc))
  ncnb = LEN_TRIM(line)
  IF (ncnb ==  0)RETURN

END IF

! strip leading blanks

done = .false.
  
DO WHILE (.NOT. done)  
  IF (line(1:1) /= ' ')THEN
    done =.TRUE.
    ncnb = ncc
  ELSE
    line(1:ncc-1)= line (2:ncc)
    ncc = ncc-1
    ncnb = ncc
  END IF
END DO

! strip double blanks

done = .false.
DO WHILE (.NOT. done)
  i = INDEX(line(1:ncc), '  ')
  IF (i == 0)THEN
    done =.TRUE.
    ncnb = ncc
  ELSE 
    line(i:ncc-1)= line (i+1:ncc)
    ncc = ncc-1
    ncnb = ncc
  END IF
END DO
! strip extraneous blanks

done = .false.
DO WHILE (.NOT. done)   
  i = INDEX(line(1:ncc), ', ')
  IF (i == 0)THEN
    done =.TRUE.
    ncnb = ncc
  ELSE
    line(i+1:ncc-1)= line (i+2:ncc)
    ncc = ncc-1
    ncnb = ncc  
  END IF
END DO

done = .false.
DO WHILE (.NOT. done)
  i = INDEX(line(1:ncc), ' ,')
  IF (i == 0)THEN
    done =.TRUE.
    ncnb = ncc
  ELSE
    line(i+1:ncc-1)= line (i+2:ncc)
    ncc = ncc-1
    ncnb = ncc
  END IF
END DO

done = .false.
DO WHILE (.NOT. done)
  i = INDEX(line(1:ncc), '.. ')
  IF (i == 0)THEN
    done =.TRUE.
    ncnb = ncc
  ELSE
    line(i+2:ncc-1)= line (i+3:ncc)
    ncc = ncc-1
    ncnb = ncc   
  END IF
END DO

done = .false.
DO WHILE (.NOT. done)
  i = INDEX(line(1:ncc), ' ..')
  IF (i == 0)THEN
    done =.TRUE.
    ncnb = ncc
  ELSE
    line(i:ncc-1)= line (i+1:ncc)
    ncc = ncc-1
    ncnb = ncc
  END IF
END DO

!strip trailing blanks    

done = .false.
DO WHILE (.NOT. done)
  if (line(ncnb:ncnb) == ' ')then
     ncnb = ncnb - 1
  else
     done = .true.
  end if
END DO


IF (LEN_TRIM(line) .EQ.  0)ncnb = 0

RETURN
END SUBROUTINE Strip_Blanks


!===================================
SUBROUTINE read_parse_mlscf (MLSCFUnit, Mlscf_data) 
!===================================

! This subroutine reads and parses MLSCF and does some preliminary checking on 
! tokens.

! Arguments

INTEGER, INTENT (IN):: MLSCFUnit 
TYPE (mlscf_T), INTENT (OUT) :: Mlscf_data

! Parameters


CHARACTER (LEN=*), PARAMETER :: lineFmt = "(a)"



! Functions


! Variables

INTEGER ::  returnStatus, perror

TYPE (mlscfKey_T), DIMENSION(MaxNoMlscfKeys) :: MlscfTable

CHARACTER (LEN=256) :: line
CHARACTER (LEN=256) :: msg
CHARACTER (LEN=2048) :: buff
CHARACTER (LEN=MlscfEntryLen) :: mnemonic, section, AnEntry
CHARACTER (LEN=maxKeyLen) :: key
CHARACTER (LEN=MaxCharValueLen)Value
INTEGER :: nc, ncnb, ncb, version, isection, i, ib, ie, igs, ios, j, jj, k
INTEGER :: CellIndex, ValueLen, keyLen, recl, status
REAL :: tmp
LOGICAL :: eof = .FALSE.
LOGICAL :: found = .FALSE.
LOGICAL :: flag_section = .FALSE.

ALLOCATE (Mlscf_data%Sections(10), STAT=status)
if (status /= 0)&
         CALL MLSMessage (MLSMSG_Error, ModuleName, 'Error allocating MLSCF')
do isection = 1, 10
   do igs = 1, MaxNoEntriesPerSection
     Mlscf_data%Sections(isection)%Entries(igs)%MlscfEntryNoKeys = 0
   end do
   Mlscf_data%Sections(isection)%NoSectionEntries = 0
end do

MlscfTable(1)%Keyword = 'InputVersionString       '
MlscfTable(1)%Type = 'string'
MlscfTable(1)%Keylen =  18
MlscfTable(2)%Keyword = 'OutputVersionString      '
MlscfTable(2)%Type = 'string'
MlscfTable(2)%Keylen =  19
MlscfTable(3)%Keyword = 'AllowClimatologyOverloads'
MlscfTable(3)%Type = 'string'
MlscfTable(3)%Keylen =  25
MlscfTable(4)%Keyword = 'source                   '
MlscfTable(4)%Type = 'string'
MlscfTable(4)%Keylen =  6
MlscfTable(5)%Keyword = 'length                   '
MlscfTable(5)%Type = 'real  '
MlscfTable(5)%Keylen =  6
MlscfTable(6)%Keyword = 'versionRange             '
MlscfTable(6)%Type = 'string'
MlscfTable(6)%Keylen =  12
MlscfTable(7)%Keyword = 'species                  '
MlscfTable(7)%Type = 'string'
MlscfTable(7)%Keylen =  7
MlscfTable(8)%Keyword = 'range                    '
MlscfTable(8)%Type = 'range '
MlscfTable(8)%Keylen =  5
MlscfTable(9)%Keyword = 'method                   '
MlscfTable(9)%Type = 'string'
MlscfTable(9)%Keylen =  6
MlscfTable(10)%Keyword = 'scale                    '
MlscfTable(10)%Type = 'real  '
MlscfTable(10)%Keylen =  5
MlscfTable(11)%Keyword ='height                   '
MlscfTable(11)%Type ='real  '
MlscfTable(11)%Keylen = 6
MlscfTable(12)%Keyword ='Overlap                  '
MlscfTable(12)%Type ='int   '
MlscfTable(12)%Keylen = 7
MlscfTable(13)%Keyword ='IdealLength              '
MlscfTable(13)%Type ='real  '
MlscfTable(13)%Keylen = 11
MlscfTable(14)%Keyword ='HomeGeodAngle            '
MlscfTable(14)%Type ='real  '
MlscfTable(14)%Keylen = 13
MlscfTable(15)%Keyword ='ScanLowerLimit           '
MlscfTable(15)%Type ='range'
MlscfTable(15)%Keylen = 14
MlscfTable(16)%Keyword ='CriticalBands            '
MlscfTable(16)%Type = 'string'
MlscfTable(16)%Keylen = 13
MlscfTable(17)%Keyword ='MaxGap                   '
MlscfTable(17)%Type ='real  '
MlscfTable(17)%Keylen = 6
MlscfTable(18)%Keyword = 'type                    '
MlscfTable(18)%Type ='string' 
MlscfTable(18)%Keylen = 4
MlscfTable(19)%Keyword ='fraction                 '
MlscfTable(19)%Type ='real  '
MlscfTable(19)%Keylen = 8
MlscfTable(20)%Keyword ='interpolationfactor      '
MlscfTable(20)%Type ='real  '
MlscfTable(20)%Keylen = 19
MlscfTable(21)%Keyword ='VersionComment           '
MlscfTable(21)%Type ='string'
MlscfTable(21)%Keylen = 14
MlscfTable(22)%Keyword ='ScanUpperLimit           '
MlscfTable(22)%Type ='range' 
MlscfTable(22)%Keylen = 14
MlscfTable(23)%Keyword ='name                     '
MlscfTable(23)%Type ='string'
MlscfTable(23)%Keylen = 4
MlscfTable(24)%Keyword ='coordinate               '
MlscfTable(24)%Type ='string'
MlscfTable(24)%Keylen = 10
MlscfTable(25)%Keyword ='values                   '
MlscfTable(25)%Type ='string'
MlscfTable(25)%Keylen = 6
MlscfTable(26)%Keyword ='vGrid                    '
MlscfTable(26)%Type ='string'
MlscfTable(26)%Keylen = 5
MlscfTable(27)%Keyword ='hGrid                    '
MlscfTable(27)%Type ='string'
MlscfTable(27)%Keylen = 5
MlscfTable(28)%Keyword ='type                     '
MlscfTable(28)%Type ='string'
MlscfTable(28)%Keylen = 4
MlscfTable(29)%Keyword ='molecule                 '
MlscfTable(29)%Type ='string'
MlscfTable(29)%Keylen = 8
MlscfTable(30)%Keyword ='unit                     '
MlscfTable(30)%Type ='string'
MlscfTable(30)%Keylen = 4
MlscfTable(31)%Keyword ='species                  '
MlscfTable(31)%Type ='string'
MlscfTable(31)%Keylen = 7
MlscfTable(32)%Keyword ='radiances                '
MlscfTable(32)%Type ='string'
MlscfTable(32)%Keylen = 9
MlscfTable(33)%Keyword ='template                 '   
MlscfTable(33)%Type ='string'
MlscfTable(33)%Keylen = 8
MlscfTable(34)%Keyword ='copy                     '
MlscfTable(34)%Type ='string'
MlscfTable(34)%Keylen = 4
MlscfTable(35)%Keyword ='l2gp_temp                '
MlscfTable(35)%Type ='string'
MlscfTable(35)%Keylen = 9
MlscfTable(36)%Keyword ='CriticalScanningModules  '
MlscfTable(36)%Type = 'string'
MlscfTable(36)%Keylen = 23
MlscfTable(37)%Keyword ='module                   '
MlscfTable(37)%Type = 'string'
MlscfTable(37)%Keylen = 6
MlscfTable(38)%Keyword ='l2aux                    '
MlscfTable(38)%Type ='string'
MlscfTable(38)%Keylen = 5

isection = 0
igs=0

! Read  input

DO WHILE (.NOT. eof)
   ncnb = 0
   DO WHILE (ios .EQ. 0 .AND. ncnb .EQ. 0)
      READ(UNIT=MLSCFUnit, IOSTAT=ios, FMT=LineFmt) line
      nc = LEN (TRIM(line))

      ! strip leading, trailing and multiple blanks from line

      CALL MLSMessage (MLSMSG_Debug, ' ', line(1:nc))

      CALL Strip_Blanks(line, nc, ncnb)

   END DO

   IF(ios .EQ. 0) THEN

      IF (INDEX(line, 'BEGIN' ) /= 0 .AND. ncnb .GT. 0) THEN
         ! Is this a BEGIN Section Line?
         isection = isection + 1
         igs = 0
         Mlscf_data%Sections(isection)%MlscfSectionName=line(7:ncnb)
         flag_section =.FALSE.
      ELSE IF (isection .EQ. 0 .AND. ncnb .GT. 0)THEN
         ! If not and there was no BEGIN before and this is not a comment line, log error and exit
         CALL MLSMessage (MLSMSG_Error, ModuleName, 'Missing BEGIN')
      ELSE

         ib = 1

         ! ignore comment lines

         IF(INDEX(line, 'END') /= 1) THEN

            buff(ib:ib+ncnb-1) = line(1:ncnb)
            ib = ib + ncnb
! copy continuation lines into buff

            DO WHILE (line(ncnb:ncnb) .EQ. '$')
               ib = ib- 1

               READ(UNIT=MLSCFUnit, IOSTAT=ios, FMT=lineFmt) line
               nc = LEN (TRIM(line))
               CALL MLSMessage (MLSMSG_Debug, ' ', line(1:nc))


               IF (ios == 0) THEN

                  CALL Strip_Blanks(line, nc, ncnb)
                  buff(ib:ib+ncnb-1) = line(1:ncnb)
                  ib = ib+ncnb
               ELSE
                  eof = .TRUE.
               END IF

            END DO
            ncb = ib-1

            ! Proces mlscf entry

            AnEntry(1:MlscfEntryLen) =' '
            ib = INDEX(buff, ',')

            IF (ib .EQ.  0) THEN
               CALL MLSMessage (MLSMSG_Error, ModuleName, "Comma expected, : "//buff(1:ncb))
            END IF

            AnEntry(1:ib-1)=buff(1:ib-1)
            ib = ib + 1
            igs = igs + 1
            Mlscf_data%Sections(isection)%Entries(igs)%MlscfEntryName = AnEntry
            CellIndex = 1

            ! Process keyword = value pairs

            DO WHILE (ib < ncb)
               ie = INDEX(buff(ib:ncb), '=')

               IF (ie .EQ. 0) THEN
                  CALL MLSMessage (MLSMSG_Error, ModuleName, "= expected, : "//buff(1:ncb))
               END IF

               ie = ie+ib-1
               key(1:maxKeyLen) = ' '
               Key = Capitalize(trim(buff(ib: ie-1))) 
               keyLen = ie-ib
               ib = ie + 1
               ie = INDEX(buff(ib:ncb), ',')

               IF (ie == 0)THEN
                  ie=ncb+1
               ELSE
                  ie=ie+ib-1
               END IF

               Value (1:MaxCharValueLen) = ' '
               Value = trim(buff(ib: ie-1))
               ValueLen = ie-ib
               ib = ie + 1

               found = .FALSE.
               i = 1

               DO WHILE ((.NOT. found) .AND. (i <= MaxNoMlscfKeys))

                  IF(key(1:KeyLen) == Capitalize(mlscfTable(i)%Keyword(1:mlscfTable(i)%Keylen)))THEN
                     found = .TRUE.
                  ELSE
                     i = i + 1
                  END IF
               END DO

               IF(found) THEN


                  Mlscf_data%Sections(isection)%Entries(igs)%MlscfEntryNoKeys = &
                       Mlscf_data%Sections(isection)%Entries(igs)%MlscfEntryNoKeys + 1
                  Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%Keyword= Key(1:Keylen)

                  IF(mlscfTable(i)%Type .EQ. 'real')THEN
                     j = 1     
                     DO WHILE((j <= ValueLen) .AND. (Value(j:j) /= ' '))
                        j = j + 1
                     END DO
                     READ (unit=Value(1:j-1), fmt=*)Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RealValue
                     Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%units=Value(j+1:ValueLen)
                  ELSE IF(mlscfTable(i)%Type .EQ. 'int')THEN
                     j = 1     
                     DO WHILE((j <= ValueLen) .AND. (Value(j:j) /= ' '))
                        j = j + 1
                     END DO
                     READ (unit=Value(1:j-1), fmt=*)Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%IntValue
                     Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%units=Value(j+1:ValueLen)
                  ELSE IF(mlscfTable(i)%Type .EQ. 'range')THEN
                     j = INDEX(Value, '..')
                     IF (j > 0)THEN
                        READ (unit=Value(1:j-1), fmt=*) &
                             Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeLowerBound
                        k=j+2

                        DO WHILE((k <= ValueLen) .AND. (Value(k:k) /= ' '))
                           k = k + 1
                        END DO

                        READ (unit=Value(j+2:k), fmt=*) &
                             Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound
                        Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%units=Value(k+1:ValueLen)

                        IF(Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound < &
                             Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound)THEN

                           tmp = Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound
                           Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeUpperBound = &
                                Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeLowerBound
                           Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%RangeLowerBound = tmp

                        END IF ! (Mlscf_data%Sections(isection)%Entries
                     ELSE
                        CALL MLSMessage (MLSMSG_Error, ModuleName, "Range expected, : "//Value)
                     END IF ! (j > 0)

                  ELSE IF(mlscfTable(i)%Type .EQ. 'string')THEN

                     j = INDEX(Value, '"')
                     IF(j > 0)THEN
                        k = INDEX(Value(j+1:ValueLen), '"')
                        IF(k > 0)THEN
                           Mlscf_data%Sections(isection)%Entries(igs)%Cells(CellIndex)%CharValue=Value(j+1:k)
                        ELSE
                           CALL MLSMessage (MLSMSG_Error, ModuleName, "Quote expected, : "//Value(1:ValueLen))
                        END IF
                     ELSE
                        CALL MLSMessage (MLSMSG_Error, ModuleName, "Quote expected, : "//Value(1:ValueLen))
                     END IF ! j > 0

                  END IF ! mlscfTable(i)%Type

               ELSE            
                  CALL MLSMessage (MLSMSG_Error, ModuleName, "Illegal Keyword : " //key(1:KeyLen))
               END IF ! (found)                      
               CellIndex = CellIndex + 1  
            END DO ! ib < ncb
         ELSE
            ! Check whether END Section matches BEGIN
            IF (Mlscf_data%Sections(isection)%MlscfSectionName .NE. line(5:ncnb))THEN
               CALL MLSMessage (MLSMSG_Warning, ModuleName, "BEGIN and END Section Names don't match " // &
                    Mlscf_data%Sections(isection)%Entries(igs)%MlscfEntryName//' , '//line(5:ncnb))
            END IF
            Mlscf_data%Sections(isection)%NoSectionEntries = igs
            Mlscf_data%NoSections= isection
            flag_section = .TRUE.
         END IF !INDEX(line, 'END') /= 1  


      END IF ! INDEX(line, 'BEGIN' ) /= 0 .and. ncnb .gt. 0

   ELSE

      eof = .TRUE.
      IF(.NOT. flag_section) THEN
         CALL MLSMessage (MLSMSG_Warning, ModuleName, 'Premature EOF &
              &          encountered')
      END IF

   END IF ! (ios .eq.  0)


END DO ! (.not. eof)



!==============================
END SUBROUTINE  read_parse_mlscf
!==============================

!=======================
END MODULE MLSCF
!=======================

!
! $Log$
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.10  2000/06/30 00:45:21  lungu
! Made MaxNoKeysPerEntry = 20.
!
