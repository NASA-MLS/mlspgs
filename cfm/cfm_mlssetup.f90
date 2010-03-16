module CFM_MLSSetup_m

   use CFM_Tree_Walker, only : Walk_Tree
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
   use MLSMessageModule, only: MLSMessageConfig
   use STRING_TABLE, only: DESTROY_CHAR_TABLE, DESTROY_HASH_TABLE, &
    & DESTROY_STRING_TABLE
   use SYMBOL_TABLE, only: DESTROY_SYMBOL_TABLE
   use Chunks_m, only: MLSChunk_T
   use L1BData, only: findMaxMaf
   use MLSFiles, only: InitializeMLSFile, mls_openFile, mls_closeFile
   use MLSCommon, only: MLSFile_T
   use SDPToolkit, only: PGS_S_SUCCESS
   use Intrinsic, only: l_hdf, l_ascii
   use Hdf, only: DFACC_RDONLY
   use input ! for l1boa file

   implicit none

   private 
   public :: CFM_MLSSetup, CFM_MLSCleanup

   character(len=20), parameter :: moduleName = "CFM_MLSSetup"
   type(MLSFile_T), dimension(:), save, pointer :: filedatabase => NULL()

   contains

   subroutine CFM_MLSSetup (retVal, ForwardModelConfigDatabase, filedatabase, fakeChunk)

      integer :: Root
      integer :: First_Section
      integer :: error
      integer, intent(out) :: retVal
      type (ForwardModelConfig_T), pointer, optional :: ForwardModelConfigDatabase(:)
      type (ForwardModelConfig_T), pointer :: dummy1(:) => NULL()
      type (MLSFile_T), dimension(:), pointer, optional :: filedatabase
      type (MLSChunk_T), intent(out), optional :: fakeChunk
      type (MLSFile_T), target :: l1bfile

      !Executables
      retVal = 0
      if (present(ForwardModelConfigDatabase)) nullify (ForwardModelConfigDatabase)

      MLSMessageConfig%useToolkit = .false.

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

      call h5open_f(error)
      if (error /= 0) then
         print *, "Error in initialize hdf5 library"
         retVal = 4
         return
      end if

      if (present(ForwardModelConfigDatabase)) then
         call Walk_Tree ( Root, First_Section, ForwardModelConfigDatabase )
      else
         call Walk_Tree ( Root, First_Section, dummy1 )
         deallocate (dummy1)
      end if

      ! Create filedatabase
      if (present(filedatabase)) then
         nullify(filedatabase)
         allocate(filedatabase(1), stat=error)
         if (error /= 0) return
         error = InitializeMLSFile(l1bfile, content='l1boa', &
            name=trim(l1boa), &
            shortName='L1BOA', type=l_hdf, access=DFACC_RDONLY)
         if (error == 0) then
            call mls_openFile(L1BFile, error)
            if (error /= PGS_S_SUCCESS ) then
               error = -1
               return
            else
               filedatabase(1) = l1bfile
            end if
            if (present(fakeChunk)) then
               fakeChunk%lastMafIndex = FindMaxMAF (filedatabase, fakeChunk%firstMafIndex )
               
               !fakeChunk%lastMafIndex = fakeChunk%lastMafIndex - fakeChunk%firstMafIndex
               ! temporary just have 1 maf, the first 5 maf in the file, belong to 
               ! previous day
               fakeChunk%lastMafIndex = 6
               fakeChunk%firstMafIndex = 6
            end if
         end if
      end if

   end subroutine CFM_MLSSetup

   subroutine CFM_MLSCleanup

      integer :: error
      integer :: i

      call destroy_char_table
      call destroy_hash_table
      call destroy_string_table
      call destroy_symbol_table
      call deallocate_decl
      call deallocate_tree
 
      do i = 1, size(filedatabase)
         call mls_closefile(filedatabase(i))
      end do
     
      if (associated(filedatabase)) deallocate(filedatabase)

      call h5close_f (error)
      if (error /= 0) then
         print *, "Error in finishing up hdf5 library"
         return
      end if

   end subroutine CFM_MLSCleanup

end module
