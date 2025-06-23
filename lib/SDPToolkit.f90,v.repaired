! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!===============================================================================
module SDPToolkit               ! F90 interface to SDP Toolkit.
!===============================================================================

  implicit NONE
  public

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! This module contains include files required to use the toolkit. Note that
! the include files with numbers in the names are generated  automagically
! when you install the toolkit. They are *not* in the toolkit tarball. If
! the toolkit has not installed properly on your system, those files will
! not be present and hence you will find that this file will  not compile.

! Plus interfaces for toolkit routines.

  include 'PGS_CBP.f'
  include 'PGS_CBP_6.f'
  include 'PGS_CSC_4.f'
  include 'PGS_EPH_5.f'
  include 'PGS_GCT.f'
  include 'PGS_IO.f'
  include 'PGS_IO_1.f'
  include 'PGS_MEM_7.f'
  include 'PGS_PC.f'
  include 'PGS_PC_9.f'
  include 'PGS_SMF.f'
  include 'PGS_TD.f'
  include 'PGS_TD_3.f'
  include 'PGS_MET_13.f'
  include 'PGS_MET.f'
  include 'PGS_DEM_14.f'
  include 'PGS_DEM.f'

! Parameters used by multiple modules

!  character (len=*), parameter :: earthModel = 'WGS84'
  character (len=*), parameter :: earthModel = 'GRS-80'

  integer, parameter :: max_orbits = 16
  integer, parameter :: spacecraftId = 6666

!    real, parameter :: Deg2Rad = 1.7453293E-02
!    real, parameter :: PI = 3.141592654
!    real, parameter :: Rad2Deg = 57.295780
   
! Warn if non-zero return from PGS_MET_Remove; toolkit currently (5.2.8.2)
! implements this in C as a void fun
! and Fortran as an int fun with 0 args
! (see $PGSHOME/src/MET/tools/PGS_MET_Remove.c and 
!  $PGSHOME/src/MET/tools/PGS_METbindFORTRAN.c)
  logical, parameter :: WARNIfCantPGSMETRemove = .FALSE.  
   
! Use the SDP Toolkit for the above functions
! (else must use a substitute or bypass calls)
! (may be reset by main program)
  logical :: UseSDPToolkit = .TRUE.

  interface

    ! Coordinate conversion
    integer function PGS_CSC_ECItoECR ( N, AsciiUTC, Offsets, &
                                      & PosVelECI, PosVelECR )
      ! Convert ECI coordinates to ECR
      integer, intent(in) :: N
      character(*), intent(in) :: AsciiUTC
      double precision, intent(in) :: Offsets(n)      ! Seconds
      double precision, intent(in) :: PosVelECI(6,n)  ! meters, meters/second
      double precision, intent(out) :: PosVelECR(6,n) ! meters, meters/second
    end function PGS_CSC_ECItoECR

    integer function PGS_CSC_ECRtoECI ( N, AsciiUTC, Offsets, &
                                      & PosVelECR, PosVelECI )
      ! Convert ECR coordinates to ECI
      integer, intent(in) :: N
      character(*), intent(in) :: AsciiUTC
      double precision, intent(in) :: Offsets(n)      ! Seconds
      double precision, intent(in) :: PosVelECR(6,n)  ! meters, meters/second
      double precision, intent(out) :: PosVelECI(6,n) ! meters, meters/second
    end function PGS_CSC_ECRtoECI

    integer function PGS_CSC_GeoToECR ( Longitude, Latitude, Altitude, &
                                      & EarthEllipseTag, POSECR )
      ! Convert geodetic coordinates to ECR
      double precision, intent(in) :: Longitude, Latitude ! Radians, geodetic
      double precision, intent(in) :: Altitude            ! Meters, geodetic
      character(*), intent(in) :: EarthEllipseTag
      double precision, intent(out) :: POSECR(3)          ! Meters
    end function PGS_CSC_GeoToECR

    ! I/O
    integer function PGS_IO_Gen_OpenF ( file_logical,file_access, &
         & record_length, file_handle, file_version )
      integer, intent(in) :: file_logical
      integer, intent(in) :: file_access
      integer, intent(in) :: record_length
      integer, intent(out) :: file_handle
      integer, intent(in) :: file_version
    end function PGS_IO_Gen_OpenF

    integer function PGS_IO_Gen_CloseF ( file_handle )
      integer, intent(in) :: file_handle
    end function PGS_IO_Gen_CloseF

    integer function PGS_PC_GetReference ( file_handle, file_version, physicalfilename )
      integer, intent(in) :: file_handle
      integer, intent(inout) :: file_version
      character (LEN=*), intent(out) :: physicalfilename
    end function PGS_PC_GetReference

    integer function PGS_PC_GetConfigData ( param_id, param_val )
      integer, intent(in) :: param_id
       character (len=*), intent(out) :: param_val
    end function PGS_PC_GetConfigData

    integer function PGS_PC_GetFileSize ( pcf_id, file_version, size )
      integer, intent(in) :: pcf_id
      integer, intent(inout) :: file_version
      integer, intent(inout) :: size
    end function PGS_PC_GetFileSize

    integer function PGS_SMF_GenerateStatusReport ( msg )
      character (LEN=*) :: msg
    end function PGS_SMF_GenerateStatusReport

    subroutine PGS_SMF_GetMsg ( CODE, MNEMONIC, MSG )
      integer, intent(out) :: CODE              ! Previously stored code
      character(len=*), intent(out) :: MNEMONIC ! Previously stored mnemonic
      character(len=*), intent(out) :: MSG      ! Previously stored message
    end subroutine PGS_SMF_GetMsg

    integer function PGS_SMF_TestStatusLevel ( status )
      integer, intent(in) :: status
    end function PGS_SMF_TestStatusLevel

    integer function PGS_TD_TAItoUTC ( sectai93, asciiutc )
      double precision, intent(in) :: sectai93
      character(LEN=27), intent(out) :: asciiutc
    end function PGS_TD_TAItoUTC

    integer function PGS_TD_asciitime_atob ( asciiutc_a, asciiutc_b )
      character(len=*), intent(in) :: asciiutc_a  ! Should be <= 27 chars
      character(len=*), intent(out) :: asciiutc_b ! Should be <= 25 chars
    end function PGS_TD_asciitime_atob

    integer function PGS_TD_asciitime_btoa ( asciiutc_b, asciiutc_a )
      character(len=*), intent(in) :: asciiutc_b  ! Should be <= 25 chars
      character(len=*), intent(out) :: asciiutc_a ! Should be <= 27 chars
    end function PGS_TD_asciitime_btoa

    integer function PGS_TD_utctotai ( time, dtime )
      character(len=*), intent(in) :: time
      double precision, intent(out) :: dtime
    end function PGS_TD_utctotai

    ! Not a toolkit routine, but one that substitutes
    ! for PGS_TD_utctotai when run w/o the pcf
    ! (Should we make a separate module to hold
    ! miscellaneous f90 interfaces such as this?)
    integer function mls_utctotai ( leapsec_file, time, dtime )
      character(len=*), intent(in) :: leapsec_file
      character(len=*), intent(in) :: time
      double precision, intent(out) :: dtime
    end function mls_utctotai

! Digital Elevation Model (DEM) functions
    integer function PGS_DEM_Open ( resolutionList, numResolutions, &
      & layerList, numLayers )
      integer, dimension(2) :: resolutionList
      integer, intent(in)   :: numResolutions
      integer, dimension(2) :: layerList
      integer, intent(in)   :: numLayers
    end function PGS_DEM_Open

    integer function PGS_DEM_Close ( resolutionList, numResolutions, &
      & layerList, numLayers )
      integer, dimension(2) :: resolutionList
      integer, intent(in)   :: numResolutions
      integer, dimension(2) :: layerList
      integer, intent(in)   :: numLayers
    end function PGS_DEM_Close

    integer function PGS_DEM_GetPoint ( resolutionList, numResolutions, &
      & layer, positionCode, &
      & pntlatitude, pntLongitude, numPoints, interpolation, interpValue )
      use ISO_C_BINDING, only: C_int16_t
      integer, dimension(2) :: resolutionList
      integer, intent(in)   :: numResolutions
      integer, intent(in)   :: positionCode
      double precision, dimension(*) :: pntlatitude
      double precision, dimension(*) :: pntlongitude
      integer, intent(in)   :: numPoints
      integer, intent(in)   :: interpolation
      integer(C_int16_t), dimension(*)   :: interpValue
    end function PGS_DEM_GetPoint

    integer function PGS_DEM_GetSize ( resolution, qualityField, &
      & positionCode, latitude, longitude, numVertPix, numHorizPix, pixByte )
      integer, intent(in)   :: resolution
      integer, intent(in)   :: qualityField
      integer, intent(in)   :: positionCode
      double precision, dimension(2) :: latitude
      double precision, dimension(2) :: longitude
      integer, intent(out)  :: numVertPix
      integer, intent(out)  :: numHorizPix
      integer, intent(out)  :: pixByte
    end function PGS_DEM_GetSize

    integer function PGS_DEM_GetQualityData ( resolution, qualityField, &
      & positionCode, latitude, longitude, qualityData )
      use ISO_C_BINDING, only: C_int16_t
      integer, intent(in)   :: resolution
      integer, intent(in)   :: qualityField
      integer, intent(in)   :: positionCode
      double precision, dimension(2) :: latitude
      double precision, dimension(2) :: longitude
      integer(C_int16_t), dimension(*)   :: qualityData
    end function PGS_DEM_GetQualityData

    integer function PGS_DEM_SortModels ( resolutionList, numResolutions, &
      & layer, positionCode, latitude, longitude, completeData )
      integer, dimension(2) :: resolutionList
      integer, intent(in)   :: numResolutions
      integer, intent(in)   :: layer
      integer, intent(in)   :: positionCode
      double precision, dimension(2) :: latitude
      double precision, dimension(2) :: longitude
      integer, intent(out)  :: completeData
    end function PGS_DEM_SortModels

! Metadata functions
! In the following, groups or imd_group will be
! character strings with length PGSd_MET_GROUP_NAME_L
! which in the current toolkit version is 105
! To see this value, execute "grep -i PGSd_MET_GROUP_NAME_L $PGSINC/*.h"
!
! groups will be an array of length PGSd_MET_NUM_OF_GROUPS
! which in the current toolkit version is 20
! To see this value, execute "grep -i PGSd_MET_NUM_OF_GROUPS $PGSINC/*.h"
!
! Unfortunately, these don't work right--don't use SDPToolkit for them
! (Shouldn't you either fix them or remove them?)
! Better declare them as externals
    integer function TK_PGS_MET_Init ( file_id, groups )
      integer, intent(in) :: file_id
      character (len=*), dimension(*) :: Groups
    end function TK_PGS_MET_Init

    integer function TK_PGS_MET_Setattr_d ( imd_group, attr_name, dval )
      character (len=*) :: imd_group
      character (len=*), intent(in) :: attr_name
      double precision, intent(in) :: dval
    end function TK_PGS_MET_Setattr_d

    integer function TK_PGS_MET_Setattr_s ( imd_group, attr_name, attr_value )
      character (len=*) :: imd_group
      character (len=*), intent(in) :: attr_name
      character (len=*), intent(in) :: attr_value
    end function TK_PGS_MET_Setattr_s

    integer function TK_PGS_MET_Getsetattr_d ( imd_group, attr_name, dval_array )
      character (len=*) :: imd_group
      character (len=*), intent(in) :: attr_name
      double precision, intent(out), DIMENSION(*) :: dval_array
    end function TK_PGS_MET_Getsetattr_d

    integer function TK_PGS_MET_Setattr_i ( imd_group, attr_name, attr_value )
      character (len=*) :: imd_group
      character (len=*), intent(in) :: attr_name
      integer, intent(in) :: attr_value
    end function TK_PGS_MET_Setattr_i

    integer function TK_PGS_MET_Write ( imd_group, hdf_attr_name, sd_id )
      character (len=*) :: imd_group
      character (len=*), intent(in) :: hdf_attr_name
      integer, intent(in) :: sd_id
    end function TK_PGS_MET_Write

    integer function TK_PGS_MET_Remove( )
    end function TK_PGS_MET_Remove

  end interface


contains 

  ! Miscellany
  ! Put here anything you want a lighter-weitght
  ! non-toolkit version of
  ! ----------------------------------------------  stamp  -----
  function stamp( chars, dateFormat, timeFormat, textCode, post, showTime )
  ! stamp input chars before outputting to PRUNIT.
  use Dates_Module, only: ReformatDate, ReformatTime
  ! Args
    character(len=*), intent(in) :: chars
    character(len=*), intent(in) :: dateFormat
    character(len=*), intent(in) :: timeFormat
    character(len=*), intent(in) :: textCode
    logical, intent(in)          :: post
    logical, intent(in)          :: showTime
    character(len=len(chars)+64) :: stamp
    !
    character(len=16) :: dateString
    character(len=16) :: timeString
    ! Executable
    stamp = chars
    if ( showTime ) then
      dateString = '' ! Intel 12 and earlier didn't fill with blanks
      timeString = '' ! Intel 12 and earlier didn't fill with blanks
      call date_and_time ( date=dateString, time=timeString )
      dateString = reFormatDate(trim(dateString), toForm=dateFormat)
      timeString = reFormatTime(trim(timeString), timeFormat)
      if ( dateFormat /= ' ' ) &
      & stamp = catStrings( stamp, dateString )
      if ( timeFormat /= ' ' ) &
      & stamp = catStrings( stamp, timeString )
    end if
    if ( textCode /= ' ' ) &
        & stamp = catStrings( stamp, textCode )
  contains
    function catStrings(a, b) result(c)
      ! Catenates strings a and b with intervening space
      ! if post then a before b
      ! otherwise b before a
      character(len=*), intent(in) :: a
      character(len=*), intent(in) :: b
      character(len = (len(a)+len(b)+1) ) :: c
      if ( post ) then
        c = trim(a) // ' ' // b
      else
        c = trim(b) // ' ' // a
      end if
    end function catStrings

  end function stamp 

!====================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module SDPToolkit
!====================

!
! $Log$
! Revision 2.24  2019/01/24 18:32:48  pwagner
! Reorganized modules that print to simplify toolkit-free builds
!
! Revision 2.23  2015/03/28 01:42:40  vsnyder
! Added some interfaces.  Deleted some uninteresting comments.  Spiffed.
!
! Revision 2.22  2014/08/06 23:04:14  vsnyder
! Use kind C_Int16_t from ISO_C_Bindings instead of INTEGER*2.  Use
! dimension(*) instead of dimension(:), the latter not being interoperable.
!
! Revision 2.21  2013/09/25 00:46:27  pwagner
! Added 2 more DEM interfaces
!
! Revision 2.20  2013/09/21 00:22:21  pwagner
! Added PGS_DEM stuff
!
! Revision 2.19  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.18  2005/12/10 00:23:56  pwagner
! Added PGS_SMF_TestStatusLevel interface
!
! Revision 2.17  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.16  2004/02/05 23:30:24  pwagner
! Flipped WARNIFCANTPGSMETREMOVE to FALSE
!
! Revision 2.15  2003/03/15 00:14:15  pwagner
! Added WARNIFCANTPGSMETREMOVE and explanation
!
! Revision 2.14  2003/03/06 19:48:08  pwagner
! Added pgs_td_asciitime_ functions
!
! Revision 2.13  2003/02/20 21:24:34  pwagner
! Added params from GCT
!
! Revision 2.12  2002/10/08 00:09:13  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.11  2002/10/01 22:03:55  pwagner
! Fixed RCS Ident Block
!
! Revision 2.10  2002/10/01 20:27:11  bwknosp
! Fixed problem
!
! Revision 2.9  2002/10/01 20:08:48  bwknosp
! Added Id and RCS info
!
! Revision 2.8  2002/04/29 17:38:10  pwagner
! Added interface for mls_utctotai
!
! Revision 2.7  2002/01/09 23:52:09  pwagner
! Removed use of MLSCommon
!
! Revision 2.6  2001/05/09 23:26:35  pwagner
! Added new functions
!
! Revision 2.5  2001/05/08 23:29:33  pwagner
! Sorry-bad metadata function interfaces
!
! Revision 2.4  2001/05/07 23:22:20  pwagner
! Added metadat functions
!
! Revision 2.3  2000/09/27 15:01:08  pumphrey
! Reinstated "numbered" Toolkit include files now I understand
! where they come from. (Duh. Groan.)
!
! Revision 2.2  2000/09/26 14:18:02  pumphrey
! Changed an arg of PGS_PC_Getreference to INOUT.
!
! Revision 2.1  2000/09/19 11:16:38  pumphrey
! Removed INCLUDEs of files that are not supplied with the SDP toolkit
!
! Revision 2.0  2000/09/05 17:41:07  dcuddy
! Change revision to 2.0
!
! Revision 1.14  2000/09/02 02:00:46  vsnyder
! Cosmetic changes
!
