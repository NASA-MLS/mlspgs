!=============================================================================
module MLSStrings               ! Some low level string handling stuff
!=============================================================================

  use MLSMessageModule

  !implicit none
  private
  public:: ReadCompleteLineWithoutComments,SplitWords,GetUniqueStrings
  public::Capitalize,depunctuate,LinearSearchStringArray
!------------------------------- RCS Ident Info ------------------------------
character(len=130),private,parameter :: id = & 
   "$Id$"
character(len=*), private,parameter ::&
     ModuleName="$RCSfile$"
!-----------------------------------------------------------------------------

!
! This module contains some low level string handling stuff
!

contains

  ! This one converts a string to all upper case (taken from HCP routine
  ! of same name) (Except that HCP can spell capitalise, that is. Fnord.)

  function Capitalize(str) result (outstr)
    ! takes a-z and replaces with A-Z 
    ! leaving other chars alone
    !--------Argument--------!
    character (len=*), intent(in) :: str
    character (len=len(str)) :: outstr

    !----------Local vars----------!
    character(len=len(str))::capstr
    integer::i,icode,offset
    !----------Executable part----------!
    capstr=str
    offset=ichar("A")-ichar("a")

    do i=1,len(str)
       icode=ichar(capstr(i:i))
       if ( icode >=ichar("a") .and. icode <= ichar("z")) then
          capstr(i:i)=char(icode+offset)
       endif
    enddo

    outstr=capstr
  end function Capitalize

  function depunctuate(str) result (outstr)
    ! Function that removes punctuation and replaces with blanks
    ! Added by HCP. This depends on the native character set being 
    ! ASCII. 
    !--------Argument--------!
    character(len=*),intent(in)::str
    character(len=len(str))::outstr
    !----------Local vars----------!
    integer::i,icode
    !----------Executable part----------!
    outstr=str
    do i=1 ,len(str)
        icode=ichar(str(i:i))
        if(  (icode >= 33 .and. icode <= 47).or.&
             (icode >= 58 .and. icode <=64).or. &
             (icode >= 91 .and. icode <=96).or. &
             (icode >= 123)) then
           outstr(i:i)=" "
        endif
    enddo

  end function depunctuate

  !---------------------------------------------------------------------------

  ! This funtion reads a line or set of lines from a text file and returns a
  ! string giving the full command with continuation lines joined, and comments
  ! removed.

  ! EOF can be returned if requested

  ! Note that this doesn't consider quotation marks, comments within quoted
  ! strings are considered comments, and continuation marks can apply within
  ! quoted strings.  Later versions of this routine may be more intelligent.

  subroutine ReadCompleteLineWithoutComments(unit,fullLine,eof, &
        commentChar,continuationChar)

    ! Dummy arguments

    integer, intent(in) :: unit ! Input file unit
    ! fullLine changed to intent InOut by HCP. Some (but not all) 
    ! F90 compilers won't let this be intent(out) because the declaration
    ! of inputLine makes use of the length of fullLine even if its _contents_
    ! are immaterial
    character(len=*), intent(inout) :: fullLine ! Output line
    character(len=*), intent(in),optional :: commentChar
    character(len=*), intent(in),optional :: continuationChar
    logical, intent(out), optional :: eof ! Set if got to eof

    ! Local variables

    integer :: ioInfo           ! IOSTAT result
    character(len=len(fullLine)) :: inputLine ! One line from file
    integer :: commentStart     ! Start of a comment in line
    integer :: lastChar         ! Last character position in line
    integer :: gotContinuation  ! 1 if continuation needed, 0 if not
    logical :: firstLine        ! A correction to be applied

    character(len=1) :: useCommentChar
    character(len=1) :: useContinuationChar

    ! Executable code

    ! Set default values for arguments
    
    if (.not. present(commentChar)) then
       useCommentChar=";"
    else
       useCommentChar=commentChar
    endif

    if (.not. present(continuationChar)) then
       useContinuationChar="$"
    else
       useContinuationChar=continuationChar
    endif

    ! Set up for loop

    fullLine=""
    firstLine=.true.
    if (present(eof))then
       eof=.false.
    endif
    readLoop: do

       ! Try to read a line of text

       read (unit=unit,fmt="(a)",iostat=ioInfo) inputLine
       if (ioInfo /= 0) then 
          if (present(eof)) then
             eof=.true.
          endif
          exit readLoop
       endif

       ! Now we look for the start of a comment and remove any following text
       ! from this line.

       commentStart=index(inputLine,useCommentChar)
       if (commentStart /= 0) then
          inputLine=inputLine(1:commentStart-1)
       endif
       ! See if the last non blank character is a contination

       lastChar=len_trim(inputLine)
       ! if bloc inserted by HCP because inputline(lastchar:lastchar:) 
       ! caused errors with array bounds checking turned on with
       ! some f90 compilers.
       if (lastChar > 0) then 
          gotContinuation=index(inputLine(lastChar:lastChar),&
               useContinuationChar)
       else
          gotContinuation=0
       endif
       ! Concatenate this with what we have so far, make sure there's an extra
       ! space there though (not for first line though)
       
       inputLine=inputLine(1:len_trim(inputLine)-gotContinuation)
       if (firstLine) then
          fullLine=inputLine
          firstLine=.false.
       else
          fullLine=fullLine(1:len_trim(fullLine)+1)//inputLine
       endif
       
       ! If we have a continuation mark, or a blank line then keep reading
       if ((gotContinuation==0).and.(len_trim(fullLine) /= 0))then
          exit readLoop
       endif
    end do readLoop

    ! Do a final trim and exit

    fullLine=trim(adjustl(fullLine))

    end subroutine ReadCompleteLineWithoutComments

    !---------------------------------------------------------------------------

    ! This subroutine is based on my IDL one of the same name.
    ! A line of input is split into sets of words.  There are two ways in which
    ! this can be invoked.  Typically it is split into `first' and 'rest'
    ! However, if the threeway option is given it is split to first, rest and
    ! last.

    ! Note that there is a slight subtlety here, spaces are treated specially
    ! due to the use of TRIM.  Thus while two commas in a row would count as
    ! two delimters, two spaces would count as one. Also if , is the delimeter
    ! then ,<space> counts as complete delimiter.

    subroutine SplitWords(line,first,rest,last,&
          threeWay,delimiter)

      ! Dummy arguments

      character (len=*), intent(in) :: line
      character (len=*), intent(out) :: first
      character (len=*), intent(out) :: rest
      character (len=*), intent(out), optional :: last

      logical, intent(in), optional :: threeWay
      character (len=*), intent(in), optional :: delimiter

      ! Local variables

      character (len=1) :: useDelimiter
      logical :: useThreeWay
      character (len=len(line))::useLine ! Line with leading spaces removed

      integer :: firstDelimiterPos,lastDelimiterPos,trimmedLen

      ! Executable code

      useLine=adjustl(line)
      trimmedLen=len_trim(useLine)

      if (present(delimiter)) then
         useDelimiter=delimiter
      else
         useDelimiter=","
      endif

      if (present(threeWay)) then
         useThreeWay=threeWay 
      else 
         useThreeWay=.false.
      endif

      ! Clear some results by default

      if (present(last))then
         last=""
      endif
      rest=""

      ! Find the first delimiter

      firstDelimiterPos=index(useLine,useDelimiter)

      if (firstDelimiterPos == 0) then
         first=useLine
      else
         first=useLine(1:firstDelimiterPos-1)
         if (useThreeWay) then
            ! In three way mode, find the last delimiter
            lastDelimiterPos=index(trim(useLine),useDelimiter,back=.true.)
            if (present(last) .and. &
                  lastDelimiterPos /= trimmedLen) then
               last=trim(useLine(lastDelimiterPos+1:))
            endif
            if (firstDelimiterPos+1 <= lastDelimiterPos-1) then
               rest=trim(useLine(firstDelimiterPos+1:lastDelimiterPos-1))
            endif
         else
            if (firstDelimiterPos /= trimmedLen) then
               rest=trim(useLine(firstDelimiterPos+1:))
            endif
         endif
      endif

     end subroutine SplitWords

     ! ------------------------------------------------------------------------
     
     ! This subroutine takes an array of strings and returns another containing
     ! only the unique entries. The resulting array is supplied by the caller
     ! Some checking is done to make sure it's appropriate
     
     subroutine GetUniqueStrings(inList,outList,noUnique)
       ! Dummy arguments
       character (len=*), dimension(:),intent(in) :: inList
       character (len=*), dimension(:),intent(out) :: outList
       integer,intent(out) :: noUnique ! Number of unique entries

       ! Local variables 
       integer :: i,j           ! Loop counters
       logical, dimension(:), allocatable :: duplicate ! Set if already found
       integer :: status        ! Status from allocate

       integer :: inSize

       ! Executable code, setup arrays

       inSize=size(inList)
       allocate (duplicate(inSize), stat=status)
       if (status /= 0) then
          i=MLSMessage(MLSMSG_Error,ModuleName, &
               MLSMSG_Allocate//"duplicate")
       endif
       do i=1,inSize
           duplicate(i)=.false.
       end do

       ! Go through and find duplicates

       do i=1,inSize-1 ! Don't bother with last one
           if (.not. duplicate(i)) then
              do j=i+1,inSize
                  if (inList(j)==inList(i)) then
                     duplicate(j)=.true.
                  endif
             end do
          endif
       end do
       
       ! Count how many unique ones there are

       noUnique=0
       do i=1,inSize
          if (.not. duplicate(i)) then
             noUnique=noUnique+1
          endif
       end do
       
       if (noUnique>size(outList)) then 
          i= MLSMessage(MLSMSG_Error,ModuleName, "outList too small")
       endif
       if (len(outList)<len(inList)) then
          i= MLSMessage(MLSMSG_Error,ModuleName, "outList strings to small")
       endif
       outList=""

       j=1
       do i=1,noUnique
          UniqueHuntLoop: do
             if (.not. duplicate(j))then
                exit UniqueHuntLoop
             endif
             j=j+1
          end do UniqueHuntLoop
          outList(i)=inList(j)
          j=j+1
       enddo

       deallocate(duplicate)
     end subroutine GetUniqueStrings

     ! ------------------------------------------------------------------------

     ! This routine does a simple linear search for a string in a list.
     ! If the case insensitive flag is set the strings are capitalized first.
     ! If the string is not found, 0 is returned

     function LinearSearchStringArray(list,string,caseInsensitive) &
          result(LinearSearchResult)
       
       ! Dummy arguments
       character (len=*),intent(in), dimension(:) :: list
       character (len=*),intent(in) :: string
       logical, optional,intent(in) :: caseInsensitive

       ! Function result
       integer :: LinearSearchResult

       ! Local variables
       integer :: i
       logical :: useCaseInsensitive
       logical :: found

       ! Executable code

       if (present(caseInsensitive)) then
          useCaseInsensitive=caseInsensitive
       else
          useCaseInsensitive=.false.
       endif

       found=.false.
       i=1

       ! Put the conditional outside the loop for speed (not that it will make 
       ! much difference for strings)

       if (useCaseInsensitive) then
          linSearchStringInsens: do
             if (Capitalize(trim(list(i))) == Capitalize(trim(string))) then
                found=.true.
                exit linSearchStringInsens
             endif
             i=i+1
             if (i>size(list)) then
                exit linSearchStringInsens
             endif
          end do linSearchStringInsens
       else
          linSearchStringSens: do
             if (trim(list(i)) == trim(string)) then
                found=.true.
                exit linSearchStringSens
             endif
             i=i+1
             if (i>size(list)) then
                exit linSearchStringSens
             endif
          end do linSearchStringSens
       endif

       if (.not. found) then
          i=0
       endif
       LinearSearchResult=i
     end function LinearSearchStringArray
       
!=============================================================================
end module MLSStrings
!=============================================================================

! $Log$
! Revision 1.14  1999/12/08 20:33:18  pumphrey
! HCP Fixed tiny bug in ReadCompleteLineWithoutComments
!
! Revision 1.13  1999/12/06 15:57:59  pumphrey
! Changed fullLine arg of ReadCompleteLineWithoutComments to InOut
!
! Revision 1.12  1999/12/06 14:19:17  pumphrey
! HCP added function depunctuate to replace non-alphanumeric chars in a
! string with blanks
!
! Revision 1.11  1999/12/01 23:01:41  livesey
! Before renaming things to upper/lower case
!
! Revision 1.10  1999/12/01 20:42:24  livesey
! Spelt RCSfile right!
!
! Revision 1.9  1999/12/01 20:42:08  livesey
! Tried File
!
! Revision 1.8  1999/12/01 20:41:04  livesey
! Tried RCSFile instead
!
! Revision 1.7  1999/12/01 20:40:22  livesey
! Changed dollar name dollar to dollar file dollar
!
! Revision 1.6  1999/12/01 20:39:18  livesey
! Changed usage of MLSMessage
!
! Revision 1.5  1999/11/11 01:38:28  livesey
! Changed indices to 1:n from 0:n-1
!
! Revision 1.4  1999/11/09 00:29:28  livesey
! Added LinearSearchStringArray
!
! Revision 1.3  1999/11/04 00:06:12  livesey
! Undid DTC's test change.
!
! Revision 1.2  1999/11/04 00:00:10  dcuddy
! dtc made change
!
! Revision 1.1  1999/11/03 23:53:57  livesey
! Added mlsstrings and signalsfile
!
