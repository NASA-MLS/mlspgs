! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module POSIX

  ! Interfaces to POSIX routines as specified by IEEE Std 1003.1-2008,
  ! 2016 edition (including corrigenda).

  ! See http://pubs.opengroup.org/onlinepubs/9699919799/

  implicit NONE
  private

  public CloseDir, FileType, GetCWD, OpenDir, ReadDir

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  interface
    function CloseDir ( DirHandle ) bind ( C, name='closedir' )
      ! Close a directory opened by OpenDir.
      use, intrinsic :: ISO_C_Binding, only: C_Ptr, C_Int
      integer(C_int) :: CloseDir
      type(C_Ptr), intent(in), value :: DirHandle
    end function CloseDir
    function FileType ( PathName, Error ) bind ( c, name='FileType' )
      ! Use this to get the file type of a symbolic link.  It needs
      ! the entire path name, with a null at the end.  If the result
      ! value is negative, the Error argument is the errno value from
      ! the C function stat(); otherwise Error is zero.
      use, intrinsic :: ISO_C_Binding, only: C_Int
      integer(C_int) :: FileType
      character(1), intent(in) :: PathName(*)
      integer(C_int), intent(out) :: Error
    end function FileType
  end interface

  ! File types provided by ReadDir and FileType, at least for Linux
  integer, parameter, public :: DT_BLK = 6      ! Block device
  integer, parameter, public :: DT_CHR = 2      ! Character device
  integer, parameter, public :: DT_DIR = 4      ! Directory
  integer, parameter, public :: DT_FIFO = 1     ! FIFO
  integer, parameter, public :: DT_LNK = 10     ! Symbolic link
  integer, parameter, public :: DT_REG = 8      ! Plain file
  integer, parameter, public :: DT_SOCK = 12    ! Socket
  integer, parameter, public :: DT_UNKNOWN = -1 ! Don't know this value

contains

  function GetCWD ( Buf, Size )
    ! Get the current working directory, and remove any C_Null_Char from it.
    use, intrinsic :: ISO_C_Binding, only: C_Associated, C_Int, C_Null_Char, &
      & C_Ptr
    type(C_Ptr) :: GetCWD ! A pointer to BUF if OK, else NULL
    character(*), intent(out) :: Buf
    integer(C_Int), value :: Size
    interface
      function POSIX_GetCWD ( Buf, Size ) bind ( C, name='getcwd' )
        use, intrinsic :: ISO_C_Binding, only: C_Ptr, C_Int
        type(C_Ptr) :: POSIX_GetCWD ! A pointer to BUF if OK, else NULL
        character(1), intent(out) :: Buf(*)
        integer(C_Int), value :: Size
      end function POSIX_GetCWD
    end interface
    integer :: I
    getCWD = POSIX_GetCWD ( Buf, Size )
    if ( c_associated(getCWD) ) then
      do
        i = index(buf,c_null_char)
        if ( i == 0 ) exit
        buf(i:) = buf(i+1:)
      end do
    end if
  end function GetCWD

  function OpenDir ( DirName )
    ! Open a directory in preparation of reading items from it using ReadDir.
    use, intrinsic :: ISO_C_Binding, only: C_Null_Char, C_Ptr
    type(C_Ptr) :: OpenDir
    character(*), intent(in) :: DirName
    interface
      function MyOpenDir ( DirName ) bind ( C, name='opendir' )
        use, intrinsic :: ISO_C_Binding, only: C_Ptr
        type(C_Ptr) :: MyOpenDir
        character(1), intent(in) :: DirName(*)
      end function MyOpenDir
    end interface
    openDir = myOpenDir ( dirName(:len_trim(dirName)) // c_null_char )
  end function OpenDir

  subroutine ReadDir ( DirName, DirHandle, NextPath, NextLen, NextType )
    ! Read next file name from the directory opened by OpenDir.  If
    ! the file type is a link, use FileType to inquire the file type to
    ! which the link refers.
    use, intrinsic :: ISO_C_Binding, only: C_Associated, C_Null_Char, &
      & C_F_Pointer, C_Ptr
    character(*), intent(in) :: DirName   ! To compose a complete path name
    type(C_Ptr), intent(in), value :: DirHandle ! Result from OpenDir
    character(*), intent(out) :: NextPath ! Complete path name of next file
    integer, intent(out) :: NextLen       ! Negative if no next file, greater
                                          ! than len(nextPath) if the complete
                                          ! path doesn't fit in NextPath, else
                                          ! len_trim(nextPath)
    integer, intent(out) :: NextType      ! File type, -1 if there is no next
                                          ! file or its complete path name
                                          ! doesn't fit in NextPath, negative
                                          ! of errno if the file type could not
                                          ! be determined, else DT_... above
    interface
      function MyReadDir ( DirHandle, FileName, FileType ) bind ( C, name='ReadDir' )
        use, intrinsic :: ISO_C_Binding, only: C_Int, C_Ptr
        integer(C_Int) :: MyReadDir ! Length of file name, -1 if there is none
        type(C_Ptr), intent(in), value :: DirHandle
        type(C_Ptr), intent(out) :: FileName    ! Pointer to file name
        integer(C_Int), intent(out) :: FileType ! Type of the file
      end function MyReadDir
    end interface
    integer :: D, Error, I, MyFileType, N
    type(C_Ptr) :: TheFile_C
    character, pointer :: TheFile_F(:)

    n = myReadDir ( dirHandle, theFile_C, myFileType )
    if ( n <= 0 .or. .not. C_associated(theFile_C) ) then
      nextPath = ''
      nextLen = -1
      nextType = -1
    else
      call c_f_pointer ( theFile_C, theFile_F, [ n ] )
      nextPath = ''
      nextLen = n
      d = len_trim(dirName)
      nextPath = dirName(:d) // "/"
      nextLen = d
      do i = 1, n
        if ( theFile_F(i) == c_null_char ) exit
        nextLen = d + 1 + i
        if ( nextLen > len(nextPath) ) then
          nextType = -1
          return
        end if
        nextPath(nextLen:nextLen) = theFile_F(i)
      end do
      nextType = myFileType
      if ( nextType == DT_Lnk ) then
        nextType = fileType ( nextPath(:nextLen) // C_Null_Char, error )
        if ( nextType < 0 ) nextType = -error
      end if
    end if
  end subroutine ReadDir

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module POSIX

! $Log$
! Revision 2.2  2016/10/21 22:58:09  vsnyder
! Make FileType interface public
!
! Revision 2.1  2016/10/21 22:53:17  vsnyder
! Initial commit
!
