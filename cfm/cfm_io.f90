! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module CFM_IO_M

   implicit none

   public :: Read_Spectroscopy, ReadDACSFilterShapes, ReadAntennaPatterns
   public :: ReadFilterShapes, ReadPointingGrids, ReadPFAFile, ReadHDF5L2PC

   private
   !---------------------------- RCS Ident Info -------------------------------
   character(len=*), parameter :: ModuleName="$RCSfile$"
   !---------------------------------------------------------------------------

   contains
   ! Read spectroscopy file and populate the spectroscopy
   ! data base.
   subroutine Read_Spectroscopy (filename, fileType)
      use MLSStrings, only: Capitalize
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error

      ! The name of the file
      character(len=*), intent(in) :: filename
      ! Only 'HDF5' is supported currently
      character(len=*), intent(in) ::filetype

      if (capitalize(fileType) == 'HDF5') then
         call Read_HDF5_Spectroscopy(filename)
      else
         call MLSMessage(MLSMSG_Error, moduleName, &
         filetype // " not supported")
      end if
   end subroutine

   subroutine Read_HDF5_Spectroscopy (filename)
       use SpectroscopyCatalog_m, only: read_spectroscopy

       character(len=*), intent(in) :: filename

       call read_spectroscopy(0, filename, 'HDF5')

   end subroutine

   ! Read DACS filter shape file and add them to the DACS filter shapes
   ! database
   subroutine ReadDACSFilterShapes (fileName)
      use FilterShapes_m, only: open_filter_shapes_file, &
                                read_DACS_filter_shapes_file, &
                                close_filter_shapes_file

      character(len=*), intent(in) :: fileName
      integer :: lun, fileIndex

      ! Executables
      call open_filter_shapes_file (fileName, lun, fileIndex)
      call read_DACS_filter_shapes_file (lun, fileIndex, 0)
      call close_filter_shapes_file ( lun )
   end subroutine

   ! Read antenna pattern file and populate the antenna pattern
   ! database
   subroutine ReadAntennaPatterns (fileName)
      use AntennaPatterns_m, only: open_antenna_patterns_file, &
                                   read_antenna_patterns_file, &
                                   close_antenna_patterns_file

      character(len=*), intent(in) :: filename
      integer :: lun

      ! Executables
      call open_antenna_patterns_file ( fileName, lun )
      call read_antenna_patterns_file ( lun, 0)
      call close_antenna_patterns_file ( lun )
   end subroutine

   ! Read filter shape file, and populate the filter shape database
   subroutine ReadFilterShapes (fileName)
      use FilterShapes_m, only: open_filter_shapes_file, &
                                read_filter_shapes_file, &
                                close_filter_shapes_file

      character(len=*), intent(in) :: filename
      integer :: lun, fileIndex

      ! Executables
      call open_filter_shapes_file ( fileName, lun, fileIndex )
      call read_filter_shapes_file ( lun, fileIndex, 0 )
      call close_filter_shapes_file ( lun )
   end subroutine

   ! Read pointing grids and populate the pointing grid database
   subroutine ReadPointingGrids (fileName)
      use PointingGrid_m, only: open_pointing_grid_file, &
                                read_pointing_grid_file, &
                                close_pointing_grid_file

      character(len=*), intent(in) :: filename
      integer :: lun

      ! Executables
      call open_pointing_grid_file ( fileName, lun )
      call read_pointing_grid_file ( lun, 0 )
      call close_pointing_grid_file ( lun )
   end subroutine

   ! Read PFA file and populate the PFA database
   subroutine ReadPFAFile (filename)
      use PFADataBase_m, only: process_PFA_File
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error

      character(len=*), intent(in) :: filename
      integer :: num

      num = process_PFA_file (filename, 0)
      if (num == 0) &
         call MLSMessage (MLSMSG_Error, moduleName, &
         "Error processing " // filename)
   end subroutine

   ! Read L2PC and populate L2PC database
   subroutine ReadHDF5L2PC (filename)
      use L2PC_m, only: ReadCompleteHDF5L2PCFile
      use MLSCommon, only: MLSFile_T
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error
      use MLSFiles, only: InitializeMLSFile, mls_openFile, mls_closeFile
      use Intrinsic, only: l_hdf
      use Hdf, only: DFACC_RDONLY

      character(len=*), intent(in) :: filename
      type(MLSFile_T), target :: file
      type(MLSFile_T), pointer :: l2pc
      integer :: error

      error = InitializeMLSFile(file, content='l2pc', &
      name=trim(filename), type=l_hdf, access=DFACC_RDONLY)
      if (error /= 0) &
         call MLSMessage (MLSMSG_Error, moduleName, &
         "Error initializing " // trim(filename))

      call mls_openFile(file, error)
      if (error /= 0) &
         call MLSMessage (MLSMSG_Error, moduleName, &
         "Error opening " // trim (filename))

      l2pc => file
      call ReadCompleteHDF5L2PCFile (l2pc, 0)

      ! The DestroyL2PCDatabase subroutine will take care
      ! of closing the file. I don't approve of this method
      ! but it's legacy code.
      !call mls_closeFile(file)
   end subroutine

   !--------------------------- end bloc --------------------------------------
   logical function not_used_here()
   character (len=*), parameter :: IdParm = &
       "$Id$"
   character (len=len(idParm)) :: Id = idParm
      not_used_here = (id(1:1) == ModuleName(1:1))
      print *, Id ! .mod files sometimes change if PRINT is added
   end function not_used_here
   !---------------------------------------------------------------------------

end module

! $Log$
! Revision 1.6  2011/11/09 17:47:44  honghanh
! Change Read_HDF5_Spectroscopy to use Read_Spectroscopy
! in SpectroscopyCatalogs_m.
!
! Revision 1.5  2010/07/08 21:39:16  honghanh
! Add ApplyBaseline to cfm_fill_m
!
! Revision 1.4  2010/06/29 16:40:23  honghanh
! Remove all function/subroutine and user type forwarding from
! all CFM modules except for from cfm.f90
!
! Revision 1.3  2010/06/29 15:53:45  honghanh
! Add copyright comments and support for CVS log in the file
!
