 ! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module MLSCommon                ! Common definitions for the MLS software
!=============================================================================

  use SDPTOOLKIT, only: PGSd_PC_FILE_PATH_MAX

  implicit none
  public

 ! === (start of toc) ===                                                 
!     c o n t e n t s                                                    
!     - - - - - - - -                                                    

!     (data types and parameters)

! i1, i2, i4    integer types
! r4, r8        floating point types
! ip, rp        integer, floating point types used in forward model
! rv            floating point type used in vector quantity values
! rm            floating point type used in matrix values
! NameLen       character-length of quantity names
! LineLen       character-length of most input
! FileNameLen   character-length of path/filenames
! BareFNLen     character-length of filenames
! L1BInfo_T     L1B data file names, etc.

!     (subroutines and functions)
! FindFirst     Find the first logical in the array that is true
! === (end of toc) ===                                                   
! === (start of api) ===
! int FindFirst (log condition(:))      
! === (end of api) ===
 private :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  character (LEN=256) :: Id = &
       "$Id$"
  character (LEN=*), parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  ! This module contains simple definitions which are common to all the MLS PGS
  ! f90 software.

  ! Firstly, these are standard numerical types, copied from HCP
  ! (again with my change in case, sorry Hugh!)

  integer, parameter:: i1=selected_int_kind(2)
  integer, parameter:: i2=selected_int_kind(4)
  integer, parameter:: i4=selected_int_kind(7)
  integer, parameter:: r4=selected_real_kind(5)
  integer, parameter:: r8=selected_real_kind(13)

  ! Now choose the precision we want by preference (may automate this through
  ! make later on, with perl or m4 or something).
  ! These are used according to the final letter in the two-letter name:
  ! Final Letter          Context           Suggested Value
  !    ----               -------           ---------------
  !     m                 Matrix            r8  (r4 to save memory)
  !     p                 Forward Model     r8
  !     v                 Vector            r8
  integer, parameter:: rm=r4
  integer, parameter:: rp=r8
  integer, parameter:: ip=i4
  integer, parameter:: rv=r8

  ! Now we have the lengths for various strings

  integer, parameter :: NameLen=32
  integer, parameter :: LineLen=132

! Shouldn't this be PGSd_PC_FILE_PATH_MAX ?
  integer, parameter :: FileNameLen=max(PGSd_PC_FILE_PATH_MAX, 132) ! was 132
  integer, parameter :: BareFNLen=64      ! Bare file name length (w/o path)

  ! --------------------------------------------------------------------------
  
  ! Moved here from MLSFiles module
  ! Information describing the files used by the mls software
  ! Stop passing file handles back & forth bewteen routines
  ! -- pass one of these instead
  type MLSFile_T
    character (LEN=8) :: type=""  ! e.g., 'ascii', 'hdf', 'swath', 'binary'
    character (LEN=8) :: access=""  ! e.g., 'rdonly', 'write', 'rdwrite'
    character (LEN=8) :: content=""  ! e.g., 'l1brad', 'l2gp', 'l2aux'
    character (LEN=FileNameLen) :: Name=""  ! its name (usu. w/path)
    integer :: File_Id=0     ! The HDF ID (handle) or io unit for the file
    integer :: PCF_Id=0      ! The PCF ID (ref), if any,  for the file
    integer :: HDFVersion=0  ! Which hdf version is the file if hdf(eos)
    logical :: StillOpen=.false.
  end type MLSFile_T

  ! The next datatype describes the information on the L1B data files in use

  type L1BInfo_T
    integer :: L1BOAId=0     ! The HDF ID (handle) for the L1BOA file
    ! Id(s) for the L1BRAD file(s)
    integer, dimension(:), pointer :: L1BRADIds=>NULL()
    character (LEN=FileNameLen) :: L1BOAFileName=""  ! L1BOA file name
    character (LEN=FileNameLen), dimension(:), pointer :: &
         & L1BRADFileNames=>NULL()
  end type L1BInfo_T

  ! --------------------------------------------------------------------------

  ! This datatype defines the `chunks' into which the input dataset is split

  type MLSChunk_T
    integer :: firstMAFIndex   ! Index of first MAF in the chunk
    integer :: lastMAFIndex    ! Index of last MAF in the chunk
    integer :: noMAFsLowerOverlap ! Number of MAFs in the lower overlap region
    integer :: noMAFsUpperOverlap ! Number of MAFs in the upper overlap region
    integer :: chunkNumber              ! Index of this chunk
    integer, dimension(:), pointer :: HGridOffsets => NULL()
    ! This for each chunk is the index of the first non-overlapped profile in 
    ! each hGrid into the relevant output (l2gp?) file.
    integer, dimension(:), pointer :: HGridTotals => NULL()
    ! This is somewhat repetetive.  It's the total number of profiles in
    ! the output hGrid.  It's only really used in parallel runs.
  end type MLSChunk_T

  ! --------------------------------------------------------------------------

  ! The TAI93 time range

  type TAI93_Range_T
    real(r8) :: startTime ! TAI93 format
    real(r8) :: endTime   ! TAI93 format
  end type TAI93_Range_T
  ! --------------------------------------------------------------------------

  contains

  ! -------------------------------------------- FindFirst --------------
  integer function FindFirst ( condition )
    ! Find the first logical in the array that is true
    logical, dimension(:), intent(in) :: CONDITION

    ! Local variables
    integer :: I                        ! Loop counter

    ! Executable code
    FindFirst = 0
    do i = 1, size(condition)
      if ( condition(i) ) then
        FindFirst = i
        return
      end if
    end do
  end function FindFirst


!=============================================================================
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSCommon
!=============================================================================

!
! $Log$
! Revision 2.17  2003/06/20 19:31:39  pwagner
! Changes to allow direct writing of products
!
! Revision 2.16  2003/02/17 03:52:49  livesey
! Bit the bullet and changed rm to r4.
!
! Revision 2.15  2002/12/05 19:44:24  pwagner
! Moved MLSFile_T from MLSFiles to MLSCommon
!
! Revision 2.14  2002/11/06 00:16:48  pwagner
! Added toc/api blocks
!
! Revision 2.13  2002/10/08 00:09:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.12  2002/09/13 18:08:12  pwagner
! May change matrix precision rm from r8
!
! Revision 2.11  2002/08/28 22:16:18  pwagner
! Added rm, rv types
!
! Revision 2.10  2002/02/19 23:10:54  pwagner
! Added BareFNLen
!
! Revision 2.9  2002/01/09 23:51:27  pwagner
! Connected FileNameLen with PGSd_PC_FILE_PATH_MAX
!
! Revision 2.8  2001/11/14 18:03:32  livesey
! Changed FindFirst to return 0 not -1 if not found
!
! Revision 2.7  2001/09/09 02:47:58  livesey
! Moved FindFirst into MLSCommon
!
! Revision 2.6.2.3  2001/09/09 01:53:27  livesey
! Bug fix
!
! Revision 2.6.2.2  2001/09/09 01:35:46  livesey
! Moved FindFirst in from MLSL2Common
!
! Revision 2.6.2.1  2001/09/08 22:32:24  livesey
! Added RP and IP
!
! Revision 2.6  2001/04/20 23:10:53  livesey
! Initialised parameters in L1BINFO
!
! Revision 2.5  2001/03/10 18:48:17  livesey
! Really nullified the pointer!
!
! Revision 2.4  2001/03/10 07:06:46  livesey
! Nullified L1BRadfileNames in L1BInfo
!
! Revision 2.3  2001/02/09 00:38:55  livesey
! Various changes
!
! Revision 2.2  2001/01/26 23:46:35  pwagner
! Restored L1BInfo from l1/MLSL1Common back to lib/MLSCommon
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.14  2000/09/02 01:58:30  vsnyder
! Cosmetic changes
!
