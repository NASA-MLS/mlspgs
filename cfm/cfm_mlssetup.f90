module CFM_MLSSetup_m

   use CFM_Tree_Walker_m, only : Walk_Tree
   use ForwardModelConfig, only: ForwardModelConfig_T, dump
   use SpectroscopyCatalog_m, only: dump
   use QuantityTemplates, only: QuantityTemplate_T
   use VGridsDatabase, only: VGrid_T
   use LEXER_CORE, only: INIT_LEXER
   use DECLARATION_TABLE, only: ALLOCATE_DECL, DEALLOCATE_DECL
   use TREE, only: ALLOCATE_TREE, DEALLOCATE_TREE
   use INIT_TABLES_MODULE, only: INIT_TABLES
   use PARSER, only: CONFIGURATION
   use TREE_CHECKER, only: CHECK_TREE
   use H5LIB, ONLY: h5open_f, h5close_f
   ! We don't message here, but don't know why I have to use
   ! MLSMessageConfig, but this is no time to find out
   use MLSMessageModule, only: MLSMessageConfig!, MLSMSG_Allocate, &
!                               MLSMessage, MLSMSG_Error, MLSMSG_Deallocate, &
!                               MLSMSG_Warning
   use STRING_TABLE, only: DESTROY_CHAR_TABLE, DESTROY_HASH_TABLE, &
    & DESTROY_STRING_TABLE
   use SYMBOL_TABLE, only: DESTROY_SYMBOL_TABLE
   use Chunks_m, only: MLSChunk_T
   use L1BData, only: findMaxMaf, NAME_LEN, READL1BDATA, L1BDATA_T, &
                      DEALLOCATEL1BDATA, AssembleL1BQtyName
   use MLSFiles, only: InitializeMLSFile, mls_openFile, &
                       mls_closeFile, GetMLSFileByType
   use MLSCommon, only: MLSFile_T, TAI93_Range_T, r8
   use SDPToolkit, only: PGS_S_SUCCESS
   use Intrinsic, only: l_hdf, l_ascii, lit_indices
   use Hdf, only: DFACC_RDONLY
   use MLSNumerics, only: Hunt
   use ChunkDivide_m, only: ChunkDivideConfig_T, &
                            GetChunkFromTimeRange => CFM_ChunkDivide
   use Init_Tables_Module, only: l_none, phyq_mafs, phyq_time, phyq_angle, &
                                 l_orbital, l_fixed, l_both, l_either, &
                                 l_ghz
   use Allocate_Deallocate, only: allocate_test, deallocate_test
   use MLSSignals_m, only: MODULES
   use String_Table, only: get_string
   use MLSFillValues, only: ISFILLVALUE
   use MLSSets, only: FINDFIRST
   use SDPToolkit, only: pgs_td_utctotai, PGS_S_SUCCESS, PGSTD_E_NO_LEAP_SECS

   implicit none

   private
   public :: CFM_MLSSetup, CFM_MLSCleanup
   public :: MLSChunk_T

   type MAFRange_T
      integer, dimension(2) :: L1BCover ! Range found in L1B files
      integer, dimension(2) :: L2Cover  ! Range covered by L2 ProcessingRange
      integer, dimension(2) :: Expanded ! L2Cover plus prior/post overlaps
   end type

   character(len=20), parameter :: moduleName = "CFM_MLSSetup"
   integer, parameter :: CCSDSLen = 27
   type(MLSFile_T), dimension(:), save, pointer :: myFiledatabase => NULL()
   type(ChunkDivideConfig_T), save :: chunkDivideConfig  ! Using default options

   contains

   ! Parse the L2CF given in standard input, and use it to read
   ! ForwardModelConfig(s), then put them in ForwardModelConfigDatabase.
   ! Also reads Spectroscopy, MLSSignals.
   ! Reads L1BOA file and put it into the filedatabase.
   ! Uses the startTime and endTime to create a MLSChunk_T object
   ! that uses by other subroutines that read data from L1BOA.
   ! Note: the use of a filedatabase, instead of just one MLSFile_T
   ! object to store L1BOA is for convenience when wrapping
   ! existing code in MLSPGS.
   ! CFM_MLSSetup is to be called only once before CFM_MLSCleanup is called.
   subroutine CFM_MLSSetup (startTime, endTime, l1boa, retVal, &
      filedatabase, fakeChunk, ForwardModelConfigDatabase)

      ! The start time of the data to be read in the format
      ! yyyy-doyThh:mm:ss.zzzz
      character(len=CCSDSlen), intent(in) :: startTime
      ! The stop time of the data to be read in the same format as startTime.
      character(len=CCSDSlen), intent(in) :: endTime
      ! The file name of L1BOA file.
      character(len=*), intent(in) :: l1boa
      ! Return value: 0 if successful, non-zero if error(s) is encountered.
      integer, intent(out) :: retVal
      ! Output: An array of MLSFile_T object, representing open file(s),
      ! including L1BOA.
      type (MLSFile_T), dimension(:), pointer :: filedatabase
      ! Output: A data holder that holds the startTime and endTime
      ! in a way that can be understood by other subroutines.
      type (MLSChunk_T), intent(out) :: fakeChunk
      ! Output: to store all ForwardModelConfig_T objects declared in L2CF.
      ! This argument will be nullified inside this subroutine.
      type (ForwardModelConfig_T), pointer, optional :: ForwardModelConfigDatabase(:)

      integer :: Root
      integer :: First_Section
      integer :: error
      type (MLSFile_T), target :: l1bfile
      type (TAI93_Range_T) :: processingRange
      type (ForwardModelConfig_T), pointer :: dummy1(:) => NULL()

      !Executables
      retVal = 0
      if (present(ForwardModelConfigDatabase)) nullify (ForwardModelConfigDatabase)

      MLSMessageConfig%useToolkit = .false.
      MLSMessageConfig%logFileUnit = -1

      call init_lexer ( n_chars=80000, n_symbols=4000, hash_table_size=611957 )
      call allocate_decl ( ndecls=8000 )
      call allocate_tree ( n_tree=2000000 )
      call init_tables
      call configuration(Root)
      if (Root <= 0) then
         print *, 'A syntax error occurred -- there is no abstract syntax tree'
         retVal = 2
         return
      end if

      call check_tree ( root, error, first_section )
      if (error /= 0) then
         print *, "Error in tree"
         retVal = 3
         return
      end if

      ! We have to call this before opening any HDF5 file
      call h5open_f(error)
      if (error /= 0) then
         print *, "Error in initialize hdf5 library"
         retVal = 4
         return
      end if

      if (present(ForwardModelConfigDatabase)) then
         nullify(ForwardModelConfigDatabase)
         call Walk_Tree ( Root, First_Section, ForwardModelConfigDatabase )
      else
         call Walk_Tree ( Root, First_Section, dummy1 )
         deallocate (dummy1)
      end if

      ! Create filedatabase
      nullify(filedatabase)
      allocate(filedatabase(1), stat=error)
      if (error /= 0) return
      myFiledatabase => filedatabase ! Remember what you allocate
                                     ! so you can deallocate later

      error = InitializeMLSFile(l1bfile, content='l1boa', &
      name=trim(l1boa), shortName='L1BOA', type=l_hdf, access=DFACC_RDONLY)
      if (error == 0) then
         call mls_openFile(L1BFile, error)
         if (error /= PGS_S_SUCCESS ) then
            error = -1
            return
         else
            filedatabase(1) = l1bfile
         end if

         ! Get the start time and end time encoded
         error = pgs_td_utctotai (startTime, processingrange%starttime)
         if (error /= PGS_S_SUCCESS .and. error /= PGSTD_E_NO_LEAP_SECS) then
            print *, "Could not convert UTC Start time to TAI"
            retVal = 1
            return
         end if

         error = pgs_td_utctotai (endtime, processingrange%endtime)
         if (error /= PGS_S_SUCCESS .and. error /= PGSTD_E_NO_LEAP_SECS) then
            print *, "Could not convert UTC End time to TAI"
            retVal = 1
            return
         end if

         ! Initialize ChunkDivideConfig
         chunkDivideConfig%method = l_orbital
         chunkDivideConfig%maxLengthFamily = phyq_mafs
         chunkDivideConfig%skipL1BCheck = .true.
         chunkDivideConfig%homeModule = l_ghz
         chunkDivideConfig%criticalModules = l_ghz
         chunkDivideConfig%homeGeodAngle = 0.0_r8
         chunkDivideConfig%maxGap = 2.0_r8
         chunkDivideConfig%maxGapFamily = phyq_mafs
         chunkDivideConfig%maxOrbY = 50000.0_r8 ! Unit is meter
         chunkDivideConfig%scanLowerLimit = (/-20000.0_r8, 10000.0_r8 /) ! unit is meter
         chunkDivideConfig%scanUpperLimit = (/ 40000.0_r8, 200000.0_r8 /) ! unit is meter
         chunkDivideConfig%lowerOverlap = 0.0
         chunkDivideConfig%upperOverlap = 0.0
         chunkDivideConfig%lowerOverlapFamily = phyq_mafs
         chunkDivideConfig%upperOverlapFamily = phyq_mafs
         chunkDivideConfig%noChunks = 1
         ! Create a fake chunk out of the start time, end time, and L1BOA
         fakeChunk = GetChunkFromTimeRange(processingRange, filedatabase, chunkDivideConfig)
      end if

   end subroutine CFM_MLSSetup

   ! Clean up allocated memory, close opened files
   subroutine CFM_MLSCleanup

      integer :: error
      integer :: i

      ! Clean up for the tree
      call destroy_char_table
      call destroy_hash_table
      call destroy_string_table
      call destroy_symbol_table
      call deallocate_decl
      call deallocate_tree

      ! Close open files
      do i = 1, size(myFiledatabase)
         call mls_closefile(myFiledatabase(i))
      end do

      ! Deallocate filedatabase
      if (associated(myFiledatabase)) deallocate(myFiledatabase)

      ! Have to call this, after we stops using HDF5 library
      call h5close_f (error)
      if (error /= 0) then
         print *, "Error in finishing up hdf5 library"
         return
      end if

   end subroutine CFM_MLSCleanup

end module
