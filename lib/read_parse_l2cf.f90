! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.


!===================================
SUBROUTINE read_l2cf (l2cf_data, returnStatus) 
!===================================

USE l2cf

IMPLICIT NONE


!------------------- RCS Ident Info -----------------------
CHARACTER(LEN=130) :: Id = &                                                    
"$Id$"
!----------------------------------------------------------

! Brief description of program
! This subroutine reads and parses L2CF and does some preliminary checking on 
! tokens.

! Arguments

TYPE (l2cf), intent = output :: L2cf_data
! Parameters


CHARACTER (LEN=*), PARAMETER :: lineFmt = "(q,<nc>a1)"

INTEGER, PARAMETER :: L2CF = 3000

! Functions

INTEGER, EXTERNAL :: Pgs_io_gen_openF, Pgs_io_gen_closeF

! Variables


CHARACTER (LEN=80) :: line
CHARACTER (LEN=500) :: msg
CHARACTER (LEN=1000) :: buff



INTEGER :: nc, returnStatus, version

LOGICAL :: eof = .FALSE.
version = 1

! Open the L2CF as a generic file for reading

returnStatus = Pgs_io_gen_openF (L2CF, PGSd_IO_Gen_RSeqFrm, 0, &
                                processL2CF, version)

IF (returnStatus /= PGS_S_SUCCESS) THEN
  call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
  print *, 'Error opening L2CF:  ', mnemonic
  print *, msg
ENDIF

! Read  input
DO WHILE (.NOT. eof)

  READ(UNIT=processL2CF, IOSTAT=ios, FMT=lineFmt) nc, line
  IF (ios /= 0) then

! strip leading, trailing and multiple blanks from line

    call STRIP_BLANKS(line, nc, ncnb)

    IF (index(line, 'BEGIN') /= 0) THEN

       section(1:) = line (7:ncnb)
    END IF

   ib = 1
   READ(UNIT=processL2CF, IOSTAT=ios, FMT=lineFmt) nc, line
   IF (ios /= 0) THEN

     call STRIP_BLANKS(line, nc, ncnb)

! ignore comment lines

     IF(index(line, ';' /= 1) THEN

       IF(index(line, 'END' /= 1) THEN

         buff(ib:ib+ncnb-1) = line(1:ncnb)
         ib = ib + ncnb
         DO WHILE (line(ncnb:ncnb) = '$')
           ib = ib- 1

           READ(UNIT=processL2CF, IOSTAT=ios, FMT=lineFmt) nc, line
           IF (ios /= 0) THEN

             call STRIP_BLANKS(line, nc, ncnb)
             buff(ib:ib+ncnb-1) = line(1:ncnb)
             ib = ib+ncnb
           ELSE
             eof = .TRUE.
           END IF

         END DO
         ncb = ib
! Proces l2cf entry
         AnEntry(1:L2cfEntryLen) =' '
         ib = index(buff, ',')
         AnEntry(1:ib-1)=buff(1:ib-1)
         ib = ib + 1
         igs = igs + 1
         L2cf_data%GlobalSettings(igs)%L2cfEntryName = AnEntry
         CellIndex = 1
! Process keyword = value pairs
         DO WHILE (ib < ncb)
           ie = index(buff(ib:ncb), '=')
           key(1:maxKeyLen) = ' '
           Key = buff(ib, ie-1) 
           keyLen = ie-ib
           ib = ie + 1
           ie = index(buff(ib:ncb), ',')
           Value (1:MaxCharValueLen) = ' '
           Value = buff(ib, ie-1)
           ValueLen = ie-ib
           ib = ie + 1
           found = .FALSE.
           i = 1
           do while (.not. found .and. i <= MaxNoL2cfKeys)
             if(key(1:KeyLen) .eq. l2cfTable(i)%Keyword(1:l2cfTable(i)%Keylen))then
               found = .true.
             else
               i = i + 1
             end if
           end do
           if(found) then
             if(section .eq. 'GlobalSettings')then
                
               L2cf_data%GlobalSettings(igs)%L2cfEntryNoKeys = &
               L2cf_data%GlobalSettings(igs)%L2cfEntryNoKeys + 1
               L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%Keyword= Key(1:Keylen)
               if(l2cfTable(i)%Type .eq. "real")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%RealValue
                 L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "int")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%IntValue
                 L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "range")then
                 j = index(Value, '..')
                 if (j > 0)then
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%RangeLowerBound
                   k=j+2
                   do while(k <= ValueLen .and. Value(k:k) /= ' ')
                    j = j + 1
                   end do              
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%RangeUpperBound
                   L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%units=Value(k:ValueLen)

                   if(L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%RangeUpperBound .lt. &
                      L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%RangeUpperBound)then

                     tmp = L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%RangeUpperBound
                     L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%RangeUpperBound = 
                     L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%RangeLowerBound
                     L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%RangeLowerBound = tmp

                   end if
                 else
                   print, *, 'Range expected, :', Value
                 end if

               else if(l2cfTable(i)%Type .eq. 'string')then

                 j = index(Value, '"')
                 if(j > 0)then
                   k = index(Value(j+1:ValueLen)
                   if(k > 0)then
                      L2cf_data%GlobalSettings(igs)%Cells(CellIndex)%CharValue=Value(j:k)
                   else
                      print, *, 'Quote expected', Value(1:ValueLen)
                   end if
                 else
                   print, *, 'Quote expected', Value(1:ValueLen)
                 end if
    
               end if
               
             else if(section .eq. 'ReadApriori')then
                
               L2cf_data%ReadApriori(ira)%L2cfEntryNoKeys = &
               L2cf_data%ReadApriori(ira)%L2cfEntryNoKeys + 1
               L2cf_data%ReadApriori(ira)%Cells(CellIndex)%Keyword= Key(1:Keylen)
               if(l2cfTable(i)%Type .eq. "real")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%ReadApriori(ira)%Cells(CellIndex)%RealValue
                 L2cf_data%ReadApriori(ira)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "int")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%ReadApriori(ira)%Cells(CellIndex)%IntValue
                 L2cf_data%ReadApriori(ira)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "range")then
                 j = index(Value, '..')
                 if (j > 0)then
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%ReadApriori(ira)%Cells(CellIndex)%RangeLowerBound
                   k=j+2
                   do while(k <= ValueLen .and. Value(k:k) /= ' ')
                    j = j + 1
                   end do              
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%ReadApriori(ira)%Cells(CellIndex)%RangeUpperBound
                   L2cf_data%ReadApriori(ira)%Cells(CellIndex)%units=Value(k:ValueLen)

                   if(L2cf_data%ReadApriori(ira)%Cells(CellIndex)%RangeUpperBound .lt. &
                      L2cf_data%ReadApriori(ira)%Cells(CellIndex)%RangeUpperBound)then

                     tmp = L2cf_data%ReadApriori(ira)%Cells(CellIndex)%RangeUpperBound
                     L2cf_data%ReadApriori(ira)%Cells(CellIndex)%RangeUpperBound = 
                     L2cf_data%ReadApriori(ira)%Cells(CellIndex)%RangeLowerBound
                     L2cf_data%ReadApriori(ira)%Cells(CellIndex)%RangeLowerBound = tmp

                   end if
                 else
                   print, *, 'Range expected, :', Value
                 end if
               else if(l2cfTable(i)%Type .eq. 'string')then

                 j = index(Value, '"')
                 if(j > 0)then
                   k = index(Value(j+1:ValueLen)
                   if(k > 0)then
                      L2cf_data%ReadApriori(ira)%Cells(CellIndex)%CharValue=Value(j:k)
                   else
                      print, *, 'Quote expected', Value(1:ValueLen)
                   end if
                 else
                   print, *, 'Quote expected', Value(1:ValueLen)
                 end if
               end if
  
             else if(section .eq. 'MergeApriori')then
                
               L2cf_data%MergeApriori(ima)%L2cfEntryNoKeys = &
               L2cf_data%MergeApriori(ima)%L2cfEntryNoKeys + 1
               L2cf_data%MergeApriori(ima)%Cells(CellIndex)%Keyword= Key(1:Keylen)
               if(l2cfTable(i)%Type .eq. "real")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%MergeApriori(ima)%Cells(CellIndex)%RealValue
                 L2cf_data%MergeApriori(ima)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "int")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%MergeApriori(ima)%Cells(CellIndex)%IntValue
                 L2cf_data%MergeApriori(ima)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "range")then
                 j = index(Value, '..')
                 if (j > 0)then
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%MergeApriori(ima)%Cells(CellIndex)%RangeLowerBound
                   k=j+2
                   do while(k <= ValueLen .and. Value(k:k) /= ' ')
                    j = j + 1
                   end do              
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%MergeApriori(ima)%Cells(CellIndex)%RangeUpperBound
                   L2cf_data%MergeApriori(ima)%Cells(CellIndex)%units=Value(k:ValueLen)

                   if(L2cf_data%MergeApriori(ima)%Cells(CellIndex)%RangeUpperBound .lt. &
                      L2cf_data%MergeApriori(ima)%Cells(CellIndex)%RangeUpperBound)then

                     tmp = L2cf_data%MergeApriori(ima)%Cells(CellIndex)%RangeUpperBound
                     L2cf_data%MergeApriori(ima)%Cells(CellIndex)%RangeUpperBound = 
                     L2cf_data%MergeApriori(ima)%Cells(CellIndex)%RangeLowerBound
                     L2cf_data%MergeApriori(ima)%Cells(CellIndex)%RangeLowerBound = tmp

                   end if
                 else
                   print, *, 'Range expected, :', Value
                 end if

               else if(l2cfTable(i)%Type .eq. 'string')then

                 j = index(Value, '"')
                 if(j > 0)then
                   k = index(Value(j+1:ValueLen), '"')
                   if(k > 0)then
                      L2cf_data%MergeAprioir(ima)%Cells(CellIndex)%CharValue=Value(j:k)
                   else
                      print, *, 'Quote expected', Value(1:ValueLen)
                   end if
                 else
                   print, *, 'Quote expected', Value(1:ValueLen)
                 end if
               end if
 
             else if(section .eq. 'ChunkDivide')then
                
               L2cf_data%ChunkDivide(icd)%L2cfEntryNoKeys = &
               L2cf_data%ChunkDivide(icd)%L2cfEntryNoKeys + 1
               L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%Keyword= Key(1:Keylen)
               if(l2cfTable(i)%Type .eq. "real")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%RealValue
                 L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "int")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%IntValue
                 L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "range")then
                 j = index(Value, '..')
                 if (j > 0)then
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%RangeLowerBound
                   k=j+2
                   do while(k <= ValueLen .and. Value(k:k) /= ' ')
                    j = j + 1
                   end do              
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%RangeUpperBound
                   L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%units=trim(Value(k:ValueLen))

                   if(L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%RangeUpperBound .lt. &
                      L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%RangeUpperBound)then

                     tmp = L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%RangeUpperBound
                     L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%RangeUpperBound = 
                     L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%RangeLowerBound
                     L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%RangeLowerBound = tmp

                   end if
                 else
                   print, *, 'Range expected, :', Value
                 end if
               else if(l2cfTable(i)%Type .eq. 'string')then

                 j = index(Value, '"')
                 if(j > 0)then
                   k = index(Value(j+1:ValueLen)
                   if(k > 0)then
                      L2cf_data%ChunkDivide(icd)%Cells(CellIndex)%CharValue=Value(j:k)
                   else
                      print, *, 'Quote expected', Value(1:ValueLen)
                   end if
                 else
                   print, *, 'Quote expected', Value(1:ValueLen)
                 end if
               end if

             else if(section .eq. 'ProfileLayout')then
                
               L2cf_data%ProfileLayout(ipl)%L2cfEntryNoKeys = &
               L2cf_data%ProfileLayout(ipl)%L2cfEntryNoKeys + 1
               L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%Keyword= trim(Key(1:Keylen))
               if(l2cfTable(i)%Type .eq. "real")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%RealValue
                 L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%units=trim(Value(j+1:ValueLen))
               else if(l2cfTable(i)%Type .eq. "int")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%IntValue
                 L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%units=trim(Value(j:ValueLen))
               else if(l2cfTable(i)%Type .eq. "range")then
                 j = index(Value, '..')
                 if (j > 0)then
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%RangeLowerBound
                   k=j+2
                   do while(k <= ValueLen .and. Value(k:k) /= ' ')
                    j = j + 1
                   end do              
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%RangeUpperBound
                   L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%units=trim(Value(k:ValueLen))

                   if(L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%RangeUpperBound .lt. &
                      L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%RangeUpperBound)then

                     tmp = L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%RangeUpperBound
                     L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%RangeUpperBound = 
                     L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%RangeLowerBound
                     L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%RangeLowerBound = tmp

                   end if
                 else
                   print, *, 'Range expected, :', Value
                 end if
               else if(l2cfTable(i)%Type .eq. 'string')then

                 j = index(Value, '"')
                 if(j > 0)then
                   k = index(Value(j+1:ValueLen)
                   if(k > 0)then
                      L2cf_data%ProfileLayout(ipl)%Cells(CellIndex)%CharValue=Value(j:k)
                   else
                      print, *, 'Quote expected', Value(1:ValueLen)
                   end if
                 else
                   print, *, 'Quote expected', Value(1:ValueLen)
                 end if

               end if

             else if(section .eq. 'Construct')then
                
               L2cf_data%Construct(ic)%L2cfEntryNoKeys = &
               L2cf_data%Construct(ic)%L2cfEntryNoKeys + 1
               L2cf_data%Construct(ic)%Cells(CellIndex)%Keyword= Key(1:Keylen)
               if(l2cfTable(i)%Type .eq. "real")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%Construct(ic)%Cells(CellIndex)%RealValue
                 L2cf_data%Construct(ic)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "int")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%Construct(ic)%Cells(CellIndex)%IntValue
                 L2cf_data%Construct(ic)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "range")then
                 j = index(Value, '..')
                 if (j > 0)then
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%Construct(ic)%Cells(CellIndex)%RangeLowerBound
                   k=j+2
                   do while(k <= ValueLen .and. Value(k:k) /= ' ')
                    j = j + 1
                   end do              
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%Construct(ic)%Cells(CellIndex)%RangeUpperBound
                   L2cf_data%Construct(ic)%Cells(CellIndex)%units=Value(k:ValueLen)

                   if(L2cf_data%Construct(ic)%Cells(CellIndex)%RangeUpperBound .lt. &
                      L2cf_data%Construct(ic)%Cells(CellIndex)%RangeUpperBound)then

                     tmp = L2cf_data%Construct(ic)%Cells(CellIndex)%RangeUpperBound
                     L2cf_data%Construct(ic)%Cells(CellIndex)%RangeUpperBound = 
                     L2cf_data%Construct(ic)%Cells(CellIndex)%RangeLowerBound
                     L2cf_data%Construct(ic)%Cells(CellIndex)%RangeLowerBound = tmp

                   end if
                 else
                   print, *, 'Range expected, :', Value
                 end if

               else if(l2cfTable(i)%Type .eq. 'string')then

                 j = index(Value, '"')
                 if(j > 0)then
                   k = index(Value(j+1:ValueLen)
                   if(k > 0)then
                      L2cf_data%Construct(ic)%Cells(CellIndex)%CharValue=Value(j:k)
                   else
                      print, *, 'Quote expected', Value(1:ValueLen)
                   end if
                 else
                   print, *, 'Quote expected', Value(1:ValueLen)
                 end if

               end if

             else if(section .eq. 'Fill')then
                
               L2cf_data%Fill(if)%L2cfEntryNoKeys = &
               L2cf_data%Fill(if)%L2cfEntryNoKeys + 1
               L2cf_data%Fill(if)%Cells(CellIndex)%Keyword= trim(Key(1:Keylen))
               if(l2cfTable(i)%Type .eq. "real")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%Fill(if)%Cells(CellIndex)%RealValue
                 L2cf_data%Fill(if)%Cells(CellIndex)%units=trim(Value(j:ValueLen))
               else if(l2cfTable(i)%Type .eq. "int")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%Fill(if)%Cells(CellIndex)%IntValue
                 L2cf_data%Fill(if)%Cells(CellIndex)%units=trim(Value(j:ValueLen))
               else if(l2cfTable(i)%Type .eq. "range")then
                 j = index(Value, '..')
                 if (j > 0)then
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%Fill(if)%Cells(CellIndex)%RangeLowerBound
                   k=j+2
                   do while(k <= ValueLen .and. Value(k:k) /= ' ')
                    j = j + 1
                   end do              
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%Fill(if)%Cells(CellIndex)%RangeUpperBound
                   L2cf_data%Fill(if)%Cells(CellIndex)%units=trim(Value(k:ValueLen))

                   if(L2cf_data%Fill(if)%Cells(CellIndex)%RangeUpperBound .lt. &
                      L2cf_data%Fill(if)%Cells(CellIndex)%RangeUpperBound)then

                     tmp = L2cf_data%Fill(if)%Cells(CellIndex)%RangeUpperBound
                     L2cf_data%Fill(if)%Cells(CellIndex)%RangeUpperBound = 
                     L2cf_data%Fill(if)%Cells(CellIndex)%RangeLowerBound
                     L2cf_data%Fill(if)%Cells(CellIndex)%RangeLowerBound = tmp

                   end if
                 else
                   print, *, 'Range expected, :', Value
                 end if
               else if(l2cfTable(i)%Type .eq. 'string')then

                 j = index(Value, '"')
                 if(j > 0)then
                   k = index(Value(j+1:ValueLen)
                   if(k > 0)then
                      L2cf_data%Fill(if)%Cells(CellIndex)%CharValue=Value(j:k)
                   else
                      print, *, 'Quote expected', Value(1:ValueLen)
                   end if
                 else
                   print, *, 'Quote expected', Value(1:ValueLen)
                 end if

               end if
           
             else if(section .eq. 'Join')then
                
               L2cf_data%Join(ij)%L2cfEntryNoKeys = &
               L2cf_data%Join(ij)%L2cfEntryNoKeys + 1
               L2cf_data%Join(ij)%Cells(CellIndex)%Keyword= Key(1:Keylen)
               if(l2cfTable(i)%Type .eq. "real")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%Join(ij)%Cells(CellIndex)%RealValue
                 L2cf_data%Join(ij)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "int")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%Join(ij)%Cells(CellIndex)%IntValue
                 L2cf_data%Join(ij)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "range")then
                 j = index(Value, '..')
                 if (j > 0)then
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%Join(ij)%Cells(CellIndex)%RangeLowerBound
                   k=j+2
                   do while(k <= ValueLen .and. Value(k:k) /= ' ')
                    j = j + 1
                   end do              
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%Join(ij)%Cells(CellIndex)%RangeUpperBound
                   L2cf_data%Join(ij)%Cells(CellIndex)%units=Value(k:ValueLen)

                   if(L2cf_data%Join(ij)%Cells(CellIndex)%RangeUpperBound .lt. &
                      L2cf_data%Join(ij)%Cells(CellIndex)%RangeUpperBound)then

                     tmp = L2cf_data%Join(ij)%Cells(CellIndex)%RangeUpperBound
                     L2cf_data%Join(ij)%Cells(CellIndex)%RangeUpperBound = 
                     L2cf_data%Join(ij)%Cells(CellIndex)%RangeLowerBound
                     L2cf_data%Join(ij)%Cells(CellIndex)%RangeLowerBound = tmp

                   end if
                 else
                   print, *, 'Range expected, :', Value
                 end if
               else if(l2cfTable(i)%Type .eq. 'string')then

                 j = index(Value, '"')
                 if(j > 0)then
                   k = index(Value(j+1:ValueLen)
                   if(k > 0)then
                      L2cf_data%join(ij)%Cells(CellIndex)%CharValue=Value(j:k)
                   else
                      print, *, 'Quote expected', Value(1:ValueLen)
                   end if
                 else
                   print, *, 'Quote expected', Value(1:ValueLen)
                 end if
               end if
              
             else if(section .eq. 'Output')then
                
               L2cf_data%Output(io)%L2cfEntryNoKeys = &
               L2cf_data%Output(io)%L2cfEntryNoKeys + 1
               L2cf_data%Output(io)%Cells(CellIndex)%Keyword= Key(1:Keylen)
               if(l2cfTable(i)%Type .eq. "real")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%Output(io)%Cells(CellIndex)%RealValue
                 L2cf_data%Output(io)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "int")then
                 j = 1     
                 do while(j <= ValueLen .and. Value(j:j) /= ' ')
                   j = j + 1
                 end do
                 read (unit=Value(1:j-1), fmt='*') L2cf_data%Output(io)%Cells(CellIndex)%IntValue
                 L2cf_data%Output(io)%Cells(CellIndex)%units=Value(j:ValueLen)
               else if(l2cfTable(i)%Type .eq. "range")then
                 j = index(Value, '..')
                 if (j > 0)then
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%Output(io)%Cells(CellIndex)%RangeLowerBound
                   k=j+2
                   do while(k <= ValueLen .and. Value(k:k) /= ' ')
                    j = j + 1
                   end do              
                   read (unit=Value(1:j-1), fmt='*') &
                   L2cf_data%Output(io)%Cells(CellIndex)%RangeUpperBound
                   L2cf_data%Output(io)%Cells(CellIndex)%units=Value(k:ValueLen)

                   if(L2cf_data%Output(io)%Cells(CellIndex)%RangeUpperBound .lt. &
                      L2cf_data%Output(io)%Cells(CellIndex)%RangeUpperBound)then

                     tmp = L2cf_data%Output(io)%Cells(CellIndex)%RangeUpperBound
                     L2cf_data%Output(io)%Cells(CellIndex)%RangeUpperBound = 
                     L2cf_data%Output(io)%Cells(CellIndex)%RangeLowerBound
                     L2cf_data%Output(io)%Cells(CellIndex)%RangeLowerBound = tmp

                   end if
                 else
                   print, *, 'Range expected, :', Value(1:ValueLen)
                 end if
               else if(l2cfTable(i)%Type .eq. 'string')then

                 j = index(Value, '"')
                 if(j > 0)then
                   k = index(Value(j+1:ValueLen)
                   if(k > 0)then
                      L2cf_data%Output(io)%Cells(CellIndex)%CharValue=Value(j:k)
                   else
                      print, *, 'Quote expected', Value(1:ValueLen)
                   end if
                 else
                   print, *, 'Quote expected', Value(1:ValueLen)
                 end if
               end if
             end if                                  
               
         END DO
       ELSE
! Check whether END Section matches BEGIN
       END IF  
     ELSE
       eof = .TRUE.
     END IF
end
end

! Close L2CF

returnStatus = Pgs_io_gen_closeF (processL2CF)

IF (returnStatus /= PGS_S_SUCCESS) THEN
  call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
  print *, 'Error closing L2CF:  ', mnemonic
  print *, msg
ENDIF

! Display the read data for verification.


! Detailed description of program
! This program uses the Toolkit to open the L2CF as a
! generic input file.  It reads the data and parses it, storing the tokens
! into the L2cf symbol table.

!===================
END PROGRAM read_l2cf
!===================

!# $Log$

