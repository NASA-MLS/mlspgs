module CFM_IO_M
   ! To be referenced from the outside
   use Read_Mie_m, only: Read_Mie
   use FilterShapes_m, only: Destroy_DACS_Filter_Database, &
                             Destroy_Filter_Shapes_Database
   use AntennaPatterns_m, only: Destroy_Ant_Patterns_Database
   use SpectroscopyCatalog_m, only: Destroy_SpectCat_Database, &
                                    Destroy_Line_Database
   use PointingGrid_m, only: Destroy_Pointing_Grid_Database
   use L2PC_m, only: DestroyL2PCDatabase
   use PFADatabase_m, only: Destroy_PFADataBase

   implicit none

   public :: Read_Spectroscopy, ReadDACSFilterShapes, ReadAntennaPatterns
   public :: ReadFilterShapes, ReadPointingGrids, ReadPFAFile, ReadHDF5L2PC
   public :: Destroy_Pointing_Grid_Database, Destroy_Ant_Patterns_Database
   public :: Destroy_DACS_Filter_Database, Destroy_Filter_Shapes_Database
   public :: Destroy_SpectCat_Database, Destroy_Line_Database
   public :: DestroyL2PCDatabase, Destroy_PFADataBase

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
      use HDF5, only: H5F_ACC_RDONLY_F, H5FOpen_F, H5FClose_F, HSize_T
      use MLSHDF5, only: LoadFromHDF5DS, LoadPtrFromHDF5DS, GetHDF5DSDims, &
                         IsHDF5DSPresent
      use SpectroscopyCatalog_m, only: Line_T, Lines, catalog_T, &
                                       MaxContinuum, MostLines, Catalog
      use MLSCommon, only: r8
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error, &
                                  MLSMSG_Allocate, MLSMSG_DeAllocate
      use MoreTree, only: GetLitIndexFromString, GetStringIndexFromString
      use Tree, only: Null_Tree
      use Intrinsic, only: L_none
      use Molecules, only: First_Molecule, Last_Molecule
      use MLSSignals_m, only: MaxSigLen
      use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
      use Parse_Signal_m, only: Parse_Signal

      character(len=*) :: filename
      ! Local variables
      integer :: iostat, fileID, Line1, nLines, lineN, j, i
      integer(hsize_t) :: Shp(1), Shp2(2) ! To get the shapes of datasets HD
      type(line_t), pointer :: MyLines(:)
      integer, pointer :: LineList(:)
      integer, pointer :: LineIndices(:)
      character(len=maxSigLen), pointer :: LineNames(:)
      integer, pointer :: QNList(:) ! Concatenation from all lines
      integer, pointer :: QNIndices(:) ! QNIndices(i) is index in
      integer, pointer :: PolarizedIndices(:) ! PolarizedIndices(i) is index in
                                   ! PolarizedList of last Polarized for line I.
      logical, pointer :: PolarizedList(:) ! Concatenation from all lines
      integer , pointer:: SignalIndices(:) ! signalIndices(i) is index in
                                 ! SidebandList and SignalList of last signal
                                 ! for line I.
      integer, pointer :: SignalList(:) ! Concatenation from all lines
      character(len=MaxSigLen) :: SignalName
      character(len=MaxSigLen), pointer :: SignalNames(:)
      integer, dimension(:), pointer :: SigInds ! From Parse_signal
      integer, pointer :: SidebandList(:) ! Concatenation from all lines
      logical :: signalError
      character(len=63) :: MoleculeName
      character(len=63), pointer :: MoleculeNames(:)
      real(r8), pointer :: Qlog(:,:)
      character(len=maxSigLen), pointer :: CatNames(:)
      real(r8), pointer :: Continuum(:,:)
      type(catalog_t), pointer :: MyCatalog(:)

      ! Executables
      signalError = .false.
      call h5fopen_f ( trim(fileName), H5F_ACC_RDONLY_F, fileID, iostat )
      if (iostat /= 0) call MLSMessage (MLSMSG_Error, ModuleName, &
         & 'Unable to open HDF5 Spectroscopy file ' // trim(fileName) // '.' )
      call getHDF5DSDims(fileID, 'Delta', shp)
      nLines = shp(1)
      line1 = 0
      if ( associated(lines) ) line1 = size(lines)
      lineN = line1 + nLines
      allocate ( myLines(lineN), stat=iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate // 'MyLines' )
      if ( associated(lines) ) then
        myLines(:line1) = lines
        deallocate ( lines, stat=iostat )
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_DeAllocate // 'Lines' )
        lines => myLines
      else
        lines => myLines
      end if
      ! Fill in the expanded part
      nullify ( lineNames, polarizedIndices, polarizedList, &
        & qnIndices, qnList, signalIndices, signalList, &
        & sidebandList, sigInds, signalNames )
      call loadPtrFromHDF5DS ( fileID, 'LineNames', lineNames )
      call loadPtrFromHDF5DS ( fileID, 'SignalNames', signalNames )

      if ( IsHDF5DSPresent ( fileID, 'PolarizedList' ) ) then
        call loadPtrFromHDF5DS ( fileID, 'PolarizedList', polarizedList )
      else
        call Allocate_test ( polarizedList, 0, 'PolarizedList', ModuleName )
      end if

      call loadPtrFromHDF5DS ( fileID, 'PolarizedIndices', polarizedIndices, &
            lowBound=line1 )
      call loadPtrFromHDF5DS ( fileID, 'QNList', qnList )
      call loadPtrFromHDF5DS ( fileID, 'QNIndices', qnIndices, lowBound=line1 )

      if ( IsHDF5DSPresent ( fileID, 'SidebandList' ) ) then
        call loadPtrFromHDF5DS ( fileID, 'SidebandList', SidebandList )
      else
        call Allocate_test ( SidebandList, 0, 'SidebandList', ModuleName )
      end if

      if ( IsHDF5DSPresent ( fileID, 'SignalList' ) ) then
        call loadPtrFromHDF5DS ( fileID, 'SignalList', SignalList )
      else
        call Allocate_test ( SignalList, 0, 'SignalList', ModuleName )
      end if

      call loadPtrFromHDF5DS ( fileID, 'SignalIndices', signalIndices, &
        & lowBound=line1 )
      call loadFromHDF5DS ( fileID, 'Delta', lines(line1+1:lineN)%delta )
      call loadFromHDF5DS ( fileID, 'EL', lines(line1+1:lineN)%el )
      call loadFromHDF5DS ( fileID, 'Gamma', lines(line1+1:lineN)%gamma )
      call loadFromHDF5DS ( fileID, 'N', lines(line1+1:lineN)%n )
      call loadFromHDF5DS ( fileID, 'N1', lines(line1+1:lineN)%n1 )
      call loadFromHDF5DS ( fileID, 'N2', lines(line1+1:lineN)%n2 )
      call loadFromHDF5DS ( fileID, 'NS', lines(line1+1:lineN)%ns )
      call loadFromHDF5DS ( fileID, 'PS', lines(line1+1:lineN)%ps )
      call loadFromHDF5DS ( fileID, 'Str', lines(line1+1:lineN)%str )
      call loadFromHDF5DS ( fileID, 'V0', lines(line1+1:lineN)%v0 )
      call loadFromHDF5DS ( fileID, 'W', lines(line1+1:lineN)%w )
      call loadFromHDF5DS ( fileID, 'UseYi', lines(line1+1:lineN)%useYi )
      do i = line1+1, lineN
        lines(i)%line_name = 0
        if ( lineNames(i) /= '' ) &
          & lines(i)%line_name = getStringIndexFromString(trim(lineNames(i)))
        ! Don't need to nullify QN, Polarized, Sidebands or Signals fields:
        ! They spring into existence nullified.
        call allocate_test ( lines(i)%qn, qnIndices(i)-qnIndices(i-1), &
          & 'Lines(i)%QN', moduleName )
        lines(i)%qn = qnList(qnIndices(i-1)+1:qnIndices(i))
        if ( signalIndices(i) /= signalIndices(i-1) ) then
          call allocate_test ( lines(i)%signals, signalIndices(i)-signalIndices(i-1), &
            & 'Lines(i)%Signals', moduleName )
          do j = 1, size(lines(i)%signals)
            call parse_signal ( trim(signalNames(signalList(signalIndices(i-1)+j))), &
              & sigInds, null_tree )
            if ( .not. associated(sigInds) ) then
              call MLSMessage (MLSMSG_error, ModuleName,  &
              & 'The string ' // trim(signalNames(signalList(signalIndices(i-1)+j))) // &
              & ' is not a signal name.' )
              signalError = .true.
            else
              lines(i)%signals(j) = sigInds(1)
            end if
            ! We can wait to deallocate sigInds until after the loop because
            ! parse_signal does allocate_test, which deallocates it first
            ! if it's allocated.
          end do ! j = 1, size(lines(i)%signals)
          call deallocate_test ( sigInds, 'SigInds', moduleName )
          call allocate_test ( lines(i)%sidebands, signalIndices(i)-signalIndices(i-1), &
            & 'Lines(i)%Sidebands', moduleName )
          lines(i)%sidebands = sidebandList(signalIndices(i-1)+1:signalIndices(i))
          if ( polarizedIndices(i) /= polarizedIndices(i-1) ) then
            call allocate_test ( lines(i)%polarized, polarizedIndices(i)-polarizedIndices(i-1), &
              & 'Lines(i)%Polarized', moduleName )
            lines(i)%polarized = polarizedList(polarizedIndices(i-1)+1:polarizedIndices(i))
          end if
        end if
      end do
      if (signalError) &
         call MLSMessage(MLSMSG_Error, moduleName, &
           'Signals in L2CF are inconsistent with signals used to create spectroscopy file')
      call deallocate_test ( lineNames, 'LineNames', moduleName )
      call deallocate_test ( signalNames, 'SignalNames', moduleName )
      call deallocate_test ( polarizedIndices, 'PolarizedIndices', moduleName )
      call deallocate_test ( polarizedList, 'PolarizedList', moduleName )
      call deallocate_test ( qnIndices, 'QNIndices', moduleName )
      call deallocate_test ( qnList, 'QNList', moduleName )
      call deallocate_test ( signalIndices, 'SignalIndices', moduleName )
      call deallocate_test ( signalList, 'SignalList', moduleName )
      call deallocate_test ( sidebandList, 'SidebandList', moduleName )
      ! Fill the catalog
      call getHDF5DSDims ( fileID, 'Continuum', shp2 )
      if ( shp2(2) /= maxContinuum ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Second dimension of continuum field of catalog has changed.' )
      allocate ( myCatalog(shp2(1)), stat=iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate // 'MyCatalog' )
      nullify ( catNames, continuum, lineIndices, lineList, moleculeNames, qlog )
      call loadPtrFromHDF5DS ( fileID, 'CatNames', catNames )
      call loadPtrFromHDF5DS ( fileID, 'Continuum', continuum )
      call loadPtrFromHDF5DS ( fileID, 'LineList', lineList )
      call loadPtrFromHDF5DS ( fileID, 'LineIndices', lineIndices, lowBound=0 )
      call loadPtrFromHDF5DS ( fileID, 'MoleculeNames', moleculeNames  )
      call loadPtrFromHDF5DS ( fileID, 'Qlog', qlog )
      call loadFromHDF5DS ( fileID, 'Mass', myCatalog%mass )
      call loadFromHDF5DS ( fileID, 'IsotopeRatio', myCatalog%defaultIsotopeRatio )
      call loadFromHDF5DS ( fileID, 'Molecule', myCatalog%molecule )
      do i = 1, size(myCatalog)
        myCatalog(i)%species_name = 0
        if ( catNames(i) /= '' ) myCatalog(i)%species_name = &
          & GetStringIndexFromString(trim(catNames(i)))
        myCatalog(i)%continuum = continuum(i,:)
        myCatalog(i)%qlog = qlog(i,:)
        if ( myCatalog(i)%molecule <= 0 ) then
          myCatalog(i)%molecule = l_none
        else
          j = getLitIndexFromString(trim(moleculeNames(myCatalog(i)%molecule)))
          if ( j < first_molecule .or. j > last_molecule ) then
            call MLSMessage ( MLSMSG_Error, moduleName, 'The string ' // &
            & trim(moleculeNames(myCatalog(i)%molecule)) // ' is not a molecule name.' )
          end if
          myCatalog(i)%molecule = j
          call allocate_test ( myCatalog(i)%lines, lineIndices(i)-lineIndices(i-1), &
            & 'MyCatalog(i)%lines', moduleName )
          myCatalog(i)%lines = lineList(lineIndices(i-1)+1:lineIndices(i)) + line1
          mostLines = max(mostLines, size(myCatalog(i)%lines))
          catalog(myCatalog(i)%molecule) = myCatalog(i)
        end if
      end do
      deallocate ( myCatalog, stat=iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Deallocate // 'MyCatalog' )
      call deallocate_test ( catNames,      'CatNames', moduleName )
      call deallocate_test ( continuum,     'Continuum', moduleName )
      call deallocate_test ( lineList,      'LineList', moduleName )
      call deallocate_test ( lineIndices,   'LineIndices', moduleName )
      call deallocate_test ( moleculeNames, 'MoleculeNames', moduleName )
      call deallocate_test ( qLog,          'QLog', moduleName )
      call H5FClose_F ( fileID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
        & 'Unable to close HDF5 Spectroscopy file ' // trim(fileName) // '.' )

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
