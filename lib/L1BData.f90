! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
module L1BData
!===============================================================================

  use Hdf, only: DFACC_READ, SFSTART
  use Lexer_Core, only: PRINT_SOURCE
  use MLSCommon, only: R8, L1BInfo_T, FileNameLen
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error, &
    & MLSMSG_L1Bread, MLSMSG_Warning
  use MoreTree, only: Get_Boolean, Get_Field_ID
  use Output_M, only: Output
  use String_Table, only: Get_String
  use TREE, only: NSONS, &
    &             SUB_ROSA, SUBTREE, DUMP_TREE_NODE, SOURCE_REF
  implicit none
  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

! Contents:

! Definition -- L1BData_T
! subroutines -- ReadL1BData
!                deallocateL1BData
!                FindL1BData

! Remarks:  This is a prototype module containing parameters, a derived type
!           definition, and a subroutine for reading the L1B data.

! Parameters

   integer, parameter :: NAME_LEN = 64		! Max len of SDS array name

! This data type is used to store quantities from an L1B data file.

  type L1BData_T
    character (LEN=NAME_LEN) :: L1BName        ! Name of field in file
    integer :: FirstMAF                        ! First major frame read
    integer :: NoMAFs                          ! # of MAFs read
    integer :: MaxMIFs                 ! Max # of MIFs/MAF in SD array
    integer :: NoAuxInds                       ! # of auxilliary indices

    integer, dimension(:), pointer :: CounterMAF => NULL() ! dimensioned (noMAFs)

    character, dimension(:,:,:), pointer :: CharField => NULL()
    real(r8), dimension(:,:,:), pointer :: DpField => NULL()
    integer, dimension(:,:,:), pointer :: IntField => NULL()
		      ! dimensioned (noAuxInds,maxMIFs,noMAFs)
  end type L1BData_T

  private ::  announce_error
  integer, private :: Error            ! Error level -- 0 = OK

contains

  subroutine L1bradSetup ( Root, L1bInfo, F_FILE, &
   & MAXNUML1BRADIDS, ILLEGALL1BRADID )
!------------------------------------------------------------------------------

! Brief description of subroutine
! Take file name from l2cf, open, and store unit no. in l1bInfo

! Arguments
    type (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset
    integer, intent(in) :: Root         ! of the l1brad file specification.
    !                                     Indexes a "spec_args" vertex.
    integer, intent(in) :: F_FILE
    integer, intent(in) :: MAXNUML1BRADIDS
    integer, intent(in) :: ILLEGALL1BRADID

    integer :: I                        ! Loop inductor, subscript
    character(len=FileNameLen) :: FileName      ! Duh
    integer :: Son                      ! Some subtree of root.
    integer :: status, sd_id
    integer, save :: ifl1=0             ! num. of L1brad files opened so far

    ! Error message codes

    error = 0

    ! Collect data from the fields. (only one legal field: file='...')

    do i = 2, nsons(root)

      son = subtree(i,root)
!      select case ( get_field_id(son) )

!      case ( f_file )
      if(get_field_id(son) == f_file) then
      
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        if(.NOT. associated(l1bInfo%L1BRADIDs)) then
          allocate ( l1bInfo%L1BRADIDs(MAXNUML1BRADIDS), stat=status )
          l1bInfo%L1BRADIDs = ILLEGALL1BRADID
          if ( status /= 0 ) &
          & call announce_error ( son, 'Allocation failed for L1BRADIDs' )
        endif
        sd_id = sfstart(Filename, DFACC_READ)
        if ( sd_id == -1 ) then
          call announce_error ( son, &
            & 'Error opening L1BRAD file: ' //Filename)
        elseif(ifl1 == MAXNUML1BRADIDS) then
          call announce_error ( son, "Cannot open any more L1BRAD files" )
          exit
        else
          ifl1 = ifl1 + 1
          l1bInfo%L1BRADIDs(ifl1) = sd_id
          l1bInfo%L1BRADFileNames(ifl1) = Filename
        end if

!      case default
      else
          call announce_error ( son, &
            & 'Unknown field specified in read l1brad' )
        ! Can't get here if the type checker worked

!      end select
      endif

    end do

!----------------------------
  end subroutine L1bradSetup
!----------------------------

!------------------------------------------------------------------------------
  subroutine L1boaSetup ( root, l1bInfo, F_FILE )
!------------------------------------------------------------------------------

! Brief description of subroutine
! Take file name from l2cf, open, and store unit no. in l1bInfo

! Arguments
    type (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset
    integer, intent(in) :: Root         ! of the l1brad file specification.
    !                                     Indexes a "spec_args" vertex.
    integer, intent(in) :: F_FILE

    integer :: I                        ! Loop inductor, subscript
    character(len=FileNameLen) :: FileName      ! Duh
    integer :: Son                      ! Some subtree of root.
    integer :: status, sd_id

    ! Error message codes

    error = 0

    ! Collect data from the fields. (only one legal field: file='...')

    do i = 2, nsons(root)

      son = subtree(i,root)
!      select case ( get_field_id(son) )

!      case ( f_file )
      if(get_field_id(son) == f_file) then

!      case ( f_file )
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        sd_id = sfstart(Filename, DFACC_READ)
        if ( sd_id == -1 ) then
          call announce_error ( son, &
            & 'Error opening L1BOA file: ' //Filename)
        else
          l1bInfo%L1BOAID = sd_id
          l1bInfo%L1BOAFileName = Filename
        end if

!      case default
      else
          call announce_error ( son, &
            & 'Unknown field specified in read l1boa' )
        ! Can't get here if the type checker worked

      endif
!      end select

    end do

!----------------------------
  end subroutine L1boaSetup
!----------------------------


!------------------------------------------------------------------------------
  subroutine ReadL1BData ( L1FileHandle, QuantityName, L1bData, NoMAFs, Flag, &
    &                       FirstMAF, LastMAF )
!------------------------------------------------------------------------------

! Brief description of subroutine

! This is a prototype for the ReadL1BData subroutine in the Open/Init task.

! Arguments

      character (len=*), intent(in) :: quantityName

      integer, intent(in) :: L1FileHandle

      integer, intent(in), optional :: firstMAF, lastMAF
 
      type( l1bdata_t ), intent(out) :: l1bData

      integer, intent(out) :: flag, noMAFs

! Parameters

      character (len=*), parameter :: INPUT_ERR = 'Error in input argument '

      integer, parameter :: MAX_VAR_DIMS = 32

! Functions

      integer, external :: sfginfo, sfn2index, sfselect, sfrdata, sfrcdata
      integer, external :: sfendacc

! Variables

      character (len=480) :: msr

      integer :: Alloc_err, Data_type, Dealloc_err, I, N_attrs, NumMAFs, Rank
      integer :: Sds1_id, Sds2_id, Sds_index, Status
      integer :: Dim_sizes(MAX_VAR_DIMS)
      integer, allocatable :: Edge(:), Start(:), Stride(:)

      logical :: FirstCheck, LastCheck

      real, allocatable :: RealField(:,:,:)

      flag = 0

! Find data sets for counterMAF & quantity by name

      sds_index = sfn2index(L1FileHandle, 'counterMAF')
      if ( sds_index == -1) call MLSMessage ( MLSMSG_Error, ModuleName, 'Failed &
                                       &to find index of counterMAF data set.')

      sds1_id = sfselect(L1FileHandle, sds_index)
      if ( sds1_id == -1) call MLSMessage ( MLSMSG_Error, ModuleName, 'Failed to &
                                     &find identifier of counterMAF data set.')

      sds_index = sfn2index(L1FileHandle, quantityName)
      if ( sds_index == -1) call MLSMessage ( MLSMSG_Error, ModuleName, 'Failed &
                                         &to find index of quantity data set.')

      sds2_id = sfselect(L1FileHandle, sds_index)
      if ( sds2_id == -1) call MLSMessage ( MLSMSG_Error, ModuleName, 'Failed to &
                             &find identifier of data set matching the index.')

! Find rank (# of dimensions), dimension sizes of quantity data set

      status = sfginfo(sds2_id, quantityName, rank, dim_sizes, data_type, &
                       n_attrs)

      if ( status == -1) call MLSMessage ( MLSMSG_Error, ModuleName, 'Failed to &
                                                      &find rank of data set.')

! allocate, based on above SD, dim info

      allocate (edge(rank), start(rank), stride(rank), STAT=alloc_err)

      if ( alloc_err /= 0 ) then
         msr = MLSMSG_allocate // '  slab dimensions.'
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if

! Set "slab" dimensions

      DO i = 1, rank
         edge(i) = dim_sizes(i)
      end do
      edge(rank) = 1

      start = 0

      stride = 1

! Fill in "indexing" values of l1b object

      l1bData%L1BName = quantityName

      if ( rank < 2 ) then
         l1bData%maxMIFs = 1
      else if ( rank > 2 ) then
         l1bData%maxMIFs = dim_sizes(2)
      else
         l1bData%maxMIFs = dim_sizes(1)
      end if

      if ( rank > 2 ) then
         l1bData%noAuxInds = dim_sizes(1)
      else
         l1bData%noAuxInds = 1
      end if

! Check input arguments, set noMAFs

      numMAFs = dim_sizes(rank)

      firstCheck = PRESENT(firstMAF)
      lastCheck = PRESENT(lastMAF)

      if ( firstCheck ) then

          if ( (firstMAF >= numMAFs) .OR. (firstMAF < 0) ) then
             msr = INPUT_ERR // 'firstMAF'
             call MLSMessage ( MLSMSG_Error, ModuleName, msr )
          end if

          l1bData%firstMAF = firstMAF

      else

          l1bData%firstMAF = 0

      end if

      if ( lastCheck ) then

          if ( lastMAF < l1bData%firstMAF ) then
             msr = INPUT_ERR // 'lastMAF'
             call MLSMessage ( MLSMSG_Error, ModuleName, msr )
          end if

          if ( lastMAF >= numMAFs ) then
             l1bData%noMAFs = numMAFs - l1bData%firstMAF
          else
             l1bData%noMAFs = lastMAF - l1bData%firstMAF + 1
          end if

      else

         l1bData%noMAFs = numMAFs - l1bData%firstMAF

      end if

      noMAFs = l1bData%noMAFs

! allocate, read counterMAF

      allocate(l1bData%counterMAF(l1bData%noMAFs), STAT=alloc_err)
      if ( alloc_err /= 0 ) then
         msr = MLSMSG_allocate // '  counterMAF pointer.'
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if

      status = sfrdata(sds1_id,  l1bData%firstMAF , stride, &
        & l1bData%noMAFs, l1bData%counterMAF )
      if ( status == -1 ) then
         msr = MLSMSG_L1BRead // 'counterMAF.'
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if

! allocate, read according to field type; nullify unused pointers

      if ( data_type == 4 ) then

! character (DFNT_CHAR8)

         allocate ( l1bData%charField(l1bData%noAuxInds,l1bData%maxMIFs,&
                   &l1bData%noMAFs), STAT=alloc_err )
         if ( alloc_err /= 0 ) then
            msr = MLSMSG_allocate // ' charField pointer.'
            call MLSMessage ( MLSMSG_Error, ModuleName, msr )
         end if
         nullify (l1bData%intField, l1bData%dpField)

         DO i = 1, l1bData%noMAFs
            start(rank) = l1bData%firstMAF + (i-1)
            status = sfrcdata( sds2_id, start, stride, edge, &
                               l1bData%charField(:,:,i) )
            if ( status == -1 ) then
               msr = MLSMSG_L1BRead // quantityName
               call MLSMessage ( MLSMSG_Error, ModuleName, msr )
            end if
         end do

      else if ( data_type == 24 ) then

! integer (DFNT_INT32)

         allocate (l1bData%intField(l1bData%noAuxInds,l1bData%maxMIFs,&
                  &l1bData%noMAFs), STAT=alloc_err)
         if ( alloc_err /= 0 ) then
            msr = MLSMSG_allocate // ' intField pointer.'
            call MLSMessage ( MLSMSG_Error, ModuleName, msr )
         end if
         nullify ( l1bData%charField, l1bData%dpField )

         DO i = 1, l1bData%noMAFs
            start(rank) = l1bData%firstMAF + (i-1)
            status = sfrdata( sds2_id, start, stride, edge, &
                              l1bData%intField(:,:,i) )
            if ( status == -1 ) then
               msr = MLSMSG_L1BRead // quantityName
               call MLSMessage ( MLSMSG_Error, ModuleName, msr )
            end if
         end do

      else

! float (REAL or DOUBLE PRECISION)

         allocate ( l1bData%dpField(l1bData%noAuxInds,l1bData%maxMIFs,&
                   & l1bData%noMAFs), STAT=alloc_err )
         if ( alloc_err /= 0 ) then
            msr = MLSMSG_allocate // ' dpField pointer.'
            call MLSMessage ( MLSMSG_Error, ModuleName, msr )
         end if
         nullify ( l1bData%charField, l1bData%intField )

! If real (DFNT_FLOAT32), allocate real variable and then move to dp pointer

         if ( data_type == 5 ) then

            allocate ( realField(l1bData%noAuxInds,l1bData%maxMIFs,&
                     & l1bData%noMAFs), STAT=alloc_err )
            if ( alloc_err /= 0 ) then
               msr = MLSMSG_allocate // ' realField pointer.'
               call MLSMessage ( MLSMSG_Error, ModuleName, msr )
            end if

            do i = 1, l1bData%noMAFs
               start(rank) = l1bData%firstMAF + (i-1)
               status = sfrdata( sds2_id, start, stride, edge, &
                                realField(:,:,i) )
               if ( status == -1 ) then
                  msr = MLSMSG_L1BRead // quantityName
                  call MLSMessage ( MLSMSG_Error, ModuleName, msr )
               end if
            end do

            l1bData%dpField = DBLE(realField)

            deallocate ( realField, STAT=dealloc_err )
            if ( dealloc_err /= 0 ) then
               call MLSMessage ( MLSMSG_Warning, ModuleName, 'Failed &
                                     &deallocation of local real variable.' )
               flag = -1
            end if


         else

! double precision (DFNT_FLOAT64)


            do i = 1, l1bData%noMAFs
               start(rank) = l1bData%firstMAF + (i-1)
               status = sfrdata(sds2_id, start, stride, edge, &
                                l1bData%dpField(:,:,i))
               if ( status == -1 ) then
                  msr = MLSMSG_L1BRead // quantityName
                  call MLSMessage ( MLSMSG_Error, ModuleName, msr )
               end if
            end do

         end if

      end if

! Terminate access to the data sets

      status = sfendacc(sds1_id)
      status = sfendacc(sds2_id)
      if ( status == -1 ) then
         call MLSMessage ( MLSMSG_Warning, ModuleName, 'Failed to terminate &
                                                     &access to data sets.' )
         flag = -1
      end if

! deallocate local variables

      deallocate (edge, start, stride, STAT=dealloc_err )

      if ( dealloc_err /= 0 ) then
         call MLSMessage ( MLSMSG_Warning, ModuleName, 'Failed deallocation of &
                                                     &local variables.' )
         flag = -1
      end if

!----------------------------
  end subroutine ReadL1BData
!----------------------------

!----------------------------------------------
  subroutine DeallocateL1BData ( l1bData, flag )
!----------------------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal field pointers of the L1BData_T
! derived type, after the calling program has finished with the data.

! Arguments

      type( L1BData_T ), intent(inout) :: L1bData

      integer, intent(out) :: Flag

! Parameters

! Functions

! Variables

      integer :: Dealloc_err

      flag = 0

      if ( associated(l1bData%counterMAF) ) deallocate ( l1bData%counterMAF, &
        &                                               STAT=dealloc_err )
      if ( associated(l1bData%charField) ) then
         deallocate ( l1bData%charField, STAT=dealloc_err )
      else if ( associated(l1bData%intField) ) then
         deallocate ( l1bData%intField, STAT=dealloc_err )
      else
         deallocate ( l1bData%dpField, STAT=dealloc_err )
      end if

      if ( dealloc_err /= 0 ) then
         call MLSMessage ( MLSMSG_Warning, ModuleName, 'Failed deallocation of &
                                                     &L1BData pointers.' )
         flag = -1
      end if

!----------------------------------
  end subroutine DeallocateL1BData
!----------------------------------

  ! ------------------------------------------- FindL1BData ----
  integer function FindL1BData ( files, fieldName )
    integer, dimension(:), intent(in) :: files ! File handles
    character (len=*), intent(in) :: fieldName ! Name of field

    ! Externals
    integer, external :: SFN2INDEX

    ! Local variables
    integer :: i

    ! Executable code
    findL1BData=0
    do i = 1, size(files)
      if ( sfn2index(files(i),fieldName) /= -1 ) then
        findL1BData = files(i)
        return
      end if
    end do
  end function FindL1BData


  ! ------------------------------------------------  announce_error  -----
  subroutine announce_error ( lcf_where, full_message, use_toolkit, &
    & error_number )
  
   ! Arguments
  
    integer, intent(in)    :: lcf_where
    character(LEN=*), intent(in)    :: full_message
    logical, intent(in), optional :: use_toolkit
    integer, intent(in), optional    :: error_number
    ! Local
!    character (len=80) :: msg, mnemonic
!    integer :: status
    logical :: just_print_it
    logical, parameter :: default_output_by_toolkit = .true.
 
    if ( present(use_toolkit) ) then
      just_print_it = .not. use_toolkit
    else if ( default_output_by_toolkit ) then
      just_print_it = .false.
    else
      just_print_it = .true.
    end if
 
    if ( .not. just_print_it ) then
      error = max(error,1)
      call output ( '***** At ' )

      if ( lcf_where > 0 ) then
          call print_source ( source_ref(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      end if

      call output ( ': ' )
      call output ( "The " );
      if ( lcf_where > 0 ) then
        call dump_tree_node ( lcf_where, 0 )
      else
        call output ( '(no lcf tree available)' )
      end if

      call output ( " Caused the following error: ", advance='yes', &
        & from_where=ModuleName )
      call output ( trim(full_message), advance='yes', &
        & from_where=ModuleName )
      if ( present(error_number) ) then
        call output ( 'error number ', advance='no' )
        call output ( error_number, places=9, advance='yes' )
      end if
    else
      print*, '***Error in module ', ModuleName
      print*, trim(full_message)
      if ( present(error_number) ) then
        print*, 'error number ', error_number
      end if
    end if

!===========================
  end subroutine announce_error
!===========================

!=================
end module L1BData
!=================

! $Log$
! Revision 2.6  2001/05/04 22:51:06  pwagner
! Now sets L1B..FileName components in ..Setup
!
! Revision 2.5  2001/05/03 23:58:38  vsnyder
! Remove conflicts from merge
!
! Revision 2.4  2001/05/03 23:16:57  livesey
! Added use statement for lexer_core
!
! Revision 2.3  2001/05/03 22:32:25  pwagner
! Added L1B..Setup for Rad and OA
!
! Revision 2.2  2001/03/03 00:06:23  livesey
! Added FindL1BData
!
! Revision 2.1  2000/10/04 01:27:34  vsnyder
! Put ONLY clauses into the USE statements
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.11  2000/06/06 21:11:55  vsnyder
! Declare EXTERNAL attribute for HDF's sf* functions
!
