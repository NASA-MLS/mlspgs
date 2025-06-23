! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module SurfaceHeight_m

  ! Read the surface height file.
  
  implicit none

  private
  ! Public procedures:
  public :: Open_Surface_Height_File
  public :: Read_Surface_Height_File
  public :: Close_Surface_Height_File

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! -----------------------------------  Open_Surface_Height_File  -----
  subroutine Open_Surface_Height_File ( TheFile )

    use MLSCommon, only: MLSFile_t
    use MLSFiles, only: MLS_OpenFile

    type(MLSFile_T), intent(inout) :: TheFile ! Created by InitializeMLSFile

    integer :: Recl, I

    inquire ( iolength=recl ) ('x', i=1, 360*6 * 180*6 * 2 ) ! 10' grid X 2 bytes
    theFile%recordLength = recl
    call MLS_OpenFile ( theFile )

  end subroutine Open_Surface_Height_File

  ! -----------------------------------  Read_Surface_Height_File  -----
  type(griddedData_t) function Read_Surface_Height_File ( TheFile, Status ) &
    & result ( theData )

    ! The surface height file is an unformatted Fortran direct-access file of
    ! 360*6 lon x180*6 lat two-byte signed two's-complement big-endian
    ! integers giving meters above mean sea level (whatever that means), with
    ! the NORTH pole and the prime meridian at the origin.  It originated from
    ! ftp://www-surf.larc.nasa.gov/pub/surf/digelev_gtopo30.map

    ! The output grid is a 180 lat x 360 lon grid giving kilometers above
    ! mean sea level with the SOUTH pole at the origin.

    use GriddedData, only: GriddedData_T, &
      & SetupNewGriddedData, NullifyGriddedData, V_is_altitude
    use Machine, only: IO_Error
    use MLSCommon, only: MLSFile_t
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error

    type(MLSFile_T), intent(in) :: TheFile   ! From Open_Surface_Height_File
    integer, intent(out), optional :: Status ! /=0 means trouble

    integer :: I, J, Stat
    integer :: imap(2160,1080)               ! 10' X 10'
    character(2) :: surfaceHeight(2160,1080) ! 2 bytes X 10' X 10'

    read ( theFile%fileId%f_id, rec=1, iostat=stat ) surfaceHeight
    if ( stat /= 0 ) then
      call io_error ( "Unable to read surface height file ", stat, trim(theFile%name) )
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to read surface height file " // trim(theFile%name) )
    else
      ! Convert two-byte big-endian signed two's-complement integers
      imap = ichar(surfaceHeight(:,:)(1:1))*256 + ichar(surfaceHeight(:,:)(2:2))
      where ( imap > 32767 ) imap = imap - 65536
      call nullifyGriddedData ( theData ) ! for Sun's still useless compiler
      ! Setup the grid
      call SetupNewGriddedData ( theData, noHeights=1, noLats=180, noLons=360 )
      theData%dateEnds = 0
      theData%dateStarts = 0
      theData%description = 'Surface height above mean sea level'
      theData%heights = 0
      theData%heightsUnits = 'km'
      theData%lsts = 0
      theData%quantityName = 'Surface height'
      theData%sourceFileName = theFile%name
      theData%szas = 0
      theData%units = 'km'
      theData%verticalCoordinate = v_is_altitude
      do i = 1, 180
        theData%lats(i) = i - 90.5
      end do
      do i = 1, 360
        theData%lons(i) = i - 0.5
      end do
      ! Put onto 1 degree X 1 degree grid using MAXVAL (don't average it).
      do i = 1, 360
        do j = 1, 180
          theData%field(1,181-j,i,1,1,1) = &
            maxval(imap((i-1)*6+1:(i-1)*6+6,(j-1)*6+1:(j-1)*6+6))
        end do ! i
      end do ! j
      where ( theData%field > 32767 ) theData%field = theData%field - 65536
      theData%field = 0.001 * theData%field ! meters -> km
    end if
    if ( present(status) ) status = stat

  end function Read_Surface_Height_File

  ! ----------------------------------  Close_Surface_Height_File  -----
  subroutine Close_Surface_Height_File ( TheFile )

    use MLSCommon, only: MLSFile_t
    use MLSFiles, only: MLS_CloseFile

    type(MLSFile_T), intent(inout) :: TheFile ! From Open_Surface_Height_File

    call MLS_CloseFile ( theFile )

  end subroutine Close_Surface_Height_File

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module SurfaceHeight_m

! $Log$
! Revision 2.4  2014/04/02 23:03:32  pwagner
! Removed redundant open_ and close_MLSFile
!
! Revision 2.3  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.2  2007/07/31 23:40:40  vsnyder
! Use constant to avoid undefined variable reference
!
! Revision 2.1  2007/01/11 20:30:10  vsnyder
! Initial commit
!
