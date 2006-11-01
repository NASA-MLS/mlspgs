! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module MLSCommon                ! Common definitions for the MLS software
!=============================================================================

  use MLSKinds ! Everything
  use SDPTOOLKIT, only: PGSd_PC_FILE_PATH_MAX

  implicit none
  private

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
! === (end of toc) ===                                                   
! === (start of api) ===
! === (end of api) ===

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains simple definitions that are common to all the MLS PGS
  ! f90 software.

  public :: FileIDs_T
  public :: InRange
  public :: MLSFile_T
  public :: L1BInfo_T
  public :: Range_T
  public :: TAI93_Range_T

  ! Make parameters gotten from MLSKinds public

  public :: i1
  public :: i2
  public :: i4
  public :: r4
  public :: r8
  public :: rm
  public :: rp
  public :: ip
  public :: rv

  ! Now we have the lengths for various strings

  integer, public, parameter :: NameLen=32
  integer, public, parameter :: LineLen=132
  integer, public, parameter :: FileNameLen=max(PGSd_PC_FILE_PATH_MAX, 132) ! was 132
  integer, public, parameter :: BareFNLen=64      ! Bare file name length (w/o path)

  real, public, parameter ::    DEFAULTUNDEFINEDVALUE = -999.99 ! Try to use in lib, l2
  ! --------------------------------------------------------------------------
  
  ! A type to hold the hdf ids

  type FileIds_T
    integer :: f_id     = 0 ! File id, handle, or io unit
    integer :: grp_id   = 0 ! group id
    integer :: sd_id    = 0 ! sd or swath id
  end type Fileids_T
  ! --------------------------------------------------------------------------

  ! A PCFid range

  type Range_T
    integer :: Bottom   = 0
    integer :: Top      = 0
  end type Range_T
  ! --------------------------------------------------------------------------

  ! Moved here from MLSFiles module
  ! Information describing the files used by the mls software
  ! Stop passing file handles back & forth between routines
  ! -- pass one of these instead
  ! (Not used yet; maybe someday)
  type MLSFile_T
    character (LEN=16) :: content=""  ! e.g., 'l1brad', 'l2gp', 'l2aux', ..
    character (LEN=8) :: lastOperation=""  ! 'open','close','read','write'
    character (LEN=FileNameLen) :: Name=""  ! its name (usu. w/path)
    character (LEN=NameLen) :: ShortName=""  ! its short name; e.g. 'H2O'
    character (LEN=8) :: typeStr=""  ! one of {'swath', 'hdf', ..}
    integer :: type=0  ! one of {l_swath, l_hdf, ..}
    integer :: access=0  ! one of {DFACC_RDONLY, DFACC_CREATE, ..}
    integer :: HDFVersion=0  ! its hdf version if hdf(eos)
    integer :: PCFId=0      ! its PCF ID (ref), if any
    integer :: recordLength=0! its max record_length, if any
    integer :: errorCode=0  ! non-zero usu. means trouble
    logical :: StillOpen=.false.
    type(Range_T) :: PCFidRange
    type(Fileids_T) :: FileID
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

  ! The TAI93 time range

  type TAI93_Range_T
    real(r8) :: startTime ! TAI93 format
    real(r8) :: endTime   ! TAI93 format
  end type TAI93_Range_T
  ! --------------------------------------------------------------------------

  contains

  elemental function inRange(arg, range) result(relation)
    ! Is arg in range?
    integer, intent(in)       :: arg
    type(Range_T), intent(in) :: range
    logical                   :: relation
    relation = (arg < (range%top + 1)) .and. (arg > (range%bottom - 1))
  end function inRange

!=============================================================================
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSCommon
!=============================================================================

!
! $Log$
! Revision 2.28  2006/11/01 20:31:38  pwagner
! House-cleaning
!
! Revision 2.27  2005/12/16 00:02:05  pwagner
! FillValue-related stuff moved to new MLSFillValues module
!
! Revision 2.26  2005/10/19 22:53:01  vsnyder
! Move kinds to MLSKinds
!
! Revision 2.25  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.24  2005/06/14 20:31:10  pwagner
! Added and changed some fields of MLSFile_T
!
! Revision 2.23  2005/05/31 17:49:15  pwagner
! Added new fields to MLSFile_T
!
! Revision 2.22  2005/05/12 20:46:46  pwagner
! Added filterValues and isFinite procedures (Should they be elsewhere?)
!
! Revision 2.21  2004/08/03 17:58:25  pwagner
! Now holds DEFAULTUNDEFINEDVALUE to be used elsewhere
!
! Revision 2.20  2004/06/10 01:00:50  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.19  2004/05/19 19:16:40  vsnyder
! Move MLSChunks_t to Chunks_m
!
! Revision 2.18  2004/01/09 00:38:04  pwagner
! Added FindNext function
!
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
