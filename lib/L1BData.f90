! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L1BData

  ! Reading and interacting with Level 1B data (HDF4)

  use dump_0, only: DUMP
  use Hdf, only: DFACC_READ, SFSTART, SFGINFO, SFN2INDEX, SFSELECT, SFRDATA_f90, &
    & SFRCDATA, SFENDACC, DFNT_CHAR8, DFNT_INT32, DFNT_FLOAT64, &
    & DFNT_FLOAT32
  use Lexer_Core, only: PRINT_SOURCE
  use MLSCommon, only: R4, R8, L1BINFO_T, FILENAMELEN
  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_ERROR, &
    & MLSMSG_L1BREAD, MLSMSG_WARNING, MLSMSG_DEALLOCATE
  use MoreTree, only: Get_Boolean, Get_Field_ID
  use Output_M, only: Output
  use String_Table, only: Get_String
  use TREE, only: NSONS, SUB_ROSA, SUBTREE, DUMP_TREE_NODE, SOURCE_REF

  implicit none

  private

  public :: L1BData_T, L1BRadSetup, L1BOASetup, ReadL1BData, DeallocateL1BData, &
    & FINDL1BDATA, NAME_LEN, PRECISIONSUFFIX, Dump

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  !---------------------------------------------------------------------------

  interface DUMP
    module procedure DumpL1BData
  end interface

  ! Parameters
  integer, parameter :: NAME_LEN = 64  ! Max len of SDS array name
  ! suffix of sd precision; check against 'grep -i precision l1/OutputL1B.f90'
  character  (len=*), parameter :: PRECISIONSUFFIX = ' precision'

  ! This data type is used to store quantities from an L1B data file.

  type L1BData_T
    character (len=name_len) :: L1BName ! Name of field in file
    integer :: FirstMAF                 ! First major frame read
    integer :: NoMAFs                   ! # of MAFs read
    integer :: MaxMIFs                  ! Max # of MIFs/MAF in SD array
    integer :: NoAuxInds                ! # of auxilliary indices

    integer, dimension(:), pointer :: CounterMAF => NULL() ! dimensioned (noMAFs)

    character, dimension(:,:,:), pointer :: CharField => NULL()
    real(r8),  dimension(:,:,:), pointer :: DpField => NULL()
    integer,   dimension(:,:,:), pointer :: IntField => NULL()
    ! all the above dimensioned (noAuxInds,maxMIFs,noMAFs)
  end type L1BData_T

  integer :: Error            ! Error level -- 0 = OK
  
  ! Error flags returned from ReadL1BData when NeverFail=TRUE
  integer, public, parameter :: NOERROR =            0
  integer, public, parameter :: NOCOUNTERMAFINDX =   NOERROR + 1
  integer, public, parameter :: NOCOUNTERMAFID =     NOCOUNTERMAFINDX + 1
  integer, public, parameter :: NOQUANTITYINDEX =    NOCOUNTERMAFID + 1
  integer, public, parameter :: NODATASETID =        NOQUANTITYINDEX + 1
  integer, public, parameter :: NODATASETRANK =      NODATASETID + 1
  integer, public, parameter :: FIRSTMAFNOTFOUND =   NODATASETRANK + 1
  integer, public, parameter :: LASTMAFNOTFOUND =    FIRSTMAFNOTFOUND + 1
  integer, public, parameter :: CANTREADCOUNTERMAF = LASTMAFNOTFOUND + 1
  integer, public, parameter :: CANTALLOCATECHARS =  CANTREADCOUNTERMAF + 1
  integer, public, parameter :: CANTREAD3DFIELD =    CANTALLOCATECHARS + 1
  integer, public, parameter :: CANTENDCOUNTERMAF =  CANTREAD3DFIELD + 1
  integer, public, parameter :: CANTENDQUANTITY =    CANTENDCOUNTERMAF + 1

contains ! ============================ MODULE PROCEDURES =======================


  ! ------------------------------------------- L1BRadSetup ------------
  subroutine L1bradSetup ( Root, L1bInfo, F_File, MaxNumL1BRadIDs, illegalL1BRadID )
    ! Take file name from l2cf, open, and store unit no. in l1bInfo
    ! Dummy arguments
    type (L1BInfo_T) :: L1BINFO         ! File handles etc. for L1B dataset
    integer, intent(in) :: ROOT         ! of the l1brad file specification.
    integer, intent(in) :: F_FILE
    integer, intent(in) :: MAXNUML1BRADIDS
    integer, intent(in) :: ILLEGALL1BRADID

    ! Local variables
    character(len=FileNameLen) :: FILENAME       ! Duh

    integer :: I                        ! Loop inductor, subscript
    integer :: SON                      ! Some subtree of root.
    integer :: STATUS                   ! Flag
    integer :: SD_ID                    ! ID from HDF

    integer, save :: IFL1 = 0           ! num. of L1brad files opened so far

    ! Exectuable code
    error = 0

    ! Collect data from the fields. (only one legal field: file='...')
    do i = 2, nsons(root)
      son = subtree(i,root)
      if(get_field_id(son) == f_file) then
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        if(.NOT. associated(l1bInfo%L1BRADIDs)) then
          allocate ( l1bInfo%L1BRADIDs(MAXNUML1BRADIDS), stat=status )
          allocate ( l1bInfo%L1BRADFileNames(MAXNUML1BRADIDS), stat=status )
          l1bInfo%L1BRADIDs = ILLEGALL1BRADID
          if ( status /= 0 ) &
            & call announce_error ( son, 'Allocation failed for l1bInfo' )
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
      else
        call announce_error ( son, &
          & 'Unknown field specified in read l1brad' )
      endif
    end do

  end subroutine L1bradSetup

  !--------------------------------------------------- L1BOASetup --------------
  subroutine L1boaSetup ( root, l1bInfo, F_FILE )
    ! Take file name from l2cf, open, and store unit no. in l1bInfo

    ! Dummy arguments
    type (L1BInfo_T) :: L1BINFO         ! File handles etc. for L1B dataset
    integer, intent(in) :: ROOT         ! of the l1brad file specification.
    integer, intent(in) :: F_FILE       ! From init_tables_module

    ! Local variables

    character(len=FileNameLen) :: FileName ! Duh

    integer :: I                        ! Loop inductor, subscript
    integer :: SON                      ! Some subtree of root.
    integer :: STATUS                   ! Flag
    integer :: SD_ID                    ! From HDF

    ! Exectuable code
    error = 0
    ! Collect data from the fields. (only one legal field: file='...')
    do i = 2, nsons(root)
      son = subtree(i,root)
      if(get_field_id(son) == f_file) then
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        sd_id = sfstart(Filename, DFACC_READ)
        if ( sd_id == -1 ) then
          call announce_error ( son, &
            & 'Error opening L1BOA file: ' //Filename)
        else
          l1bInfo%L1BOAID = sd_id
          l1bInfo%L1BOAFileName = Filename
        end if
      else
        call announce_error ( son, &
          & 'Unknown field specified in read l1boa' )
      endif
    end do
  end subroutine L1boaSetup

  !---------------------------------------------------- ReadL1BData -------------
  subroutine ReadL1BData ( L1FileHandle, QuantityName, L1bData, NoMAFs, Flag, &
    & FirstMAF, LastMAF, NEVERFAIL )
    
    ! Dummy arguments
    character(len=*), intent(in)   :: QUANTITYNAME ! Name of SD to read
    integer, intent(in)            :: L1FILEHANDLE ! From HDF
    integer, intent(in), optional  :: FIRSTMAF ! First to read (default 0)
    integer, intent(in), optional  :: LASTMAF ! Last to read (default last in file)
    logical, intent(in), optional  :: NEVERFAIL ! Don't call MLSMessage if TRUE
    type(l1bdata_t), intent(out)   :: L1BDATA ! Result
    integer, intent(out) :: FLAG        ! Error flag
    integer, intent(out) :: NOMAFS      ! Number actually read

    ! Local Parameters
    character (len=*), parameter :: INPUT_ERR = 'Error in input argument '
    integer, parameter :: MAX_VAR_DIMS = 32

    ! Local Variables

    character (len=128) :: DUMMY        ! Dummy quantity name

    integer :: ALLOC_ERR
    integer :: DATA_TYPE
    integer :: DIM_SIZES(MAX_VAR_DIMS)
    integer :: I
    logical :: MyNeverFail
    integer :: N_ATTRS
    integer :: NUMMAFS
    integer :: RANK
    integer :: SDS1_ID
    integer :: SDS2_ID
    integer :: SDS_INDEX
    integer :: STATUS

    integer, dimension(:), pointer :: EDGE
    integer, dimension(:), pointer :: START
    integer, dimension(:), pointer :: STRIDE

    real(r4), pointer, dimension(:,:,:) :: tmpR4Field

    ! Executable code
    nullify ( edge, start, stride, tmpR4Field )
    flag = 0
    MyNeverFail = .false.
    if ( present(NeverFail) ) MyNeverFail = NeverFail

    ! Find data sets for counterMAF & quantity by name

    sds_index = sfn2index(L1FileHandle, 'counterMAF')
    if ( sds_index == -1) then
      flag = NOCOUNTERMAFINDX
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Failed to find index of counterMAF data set.')
    endif

    sds1_id = sfselect(L1FileHandle, sds_index)
    if ( sds1_id == -1) then
      flag = NOCOUNTERMAFID
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Failed to find identifier of counterMAF data set.')
    endif

    sds_index = sfn2index(L1FileHandle, quantityName)
    if ( sds_index == -1) then
      flag = NOQUANTITYINDEX
      if ( MyNeverFail ) return
      dummy = 'Failed to find index of quantity "' // trim(quantityName) // &
        & '" data set.'
      call MLSMessage ( MLSMSG_Error, ModuleName, dummy )
    end if

    sds2_id = sfselect(L1FileHandle, sds_index)
    if ( sds2_id == -1) then
      flag = NODATASETID
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Failed to find identifier of data set matching the index.')
    endif

    ! Find rank (# of dimensions), dimension sizes of quantity data set
    status = sfginfo ( sds2_id, dummy, rank, dim_sizes, data_type, &
      n_attrs )

    if ( status == -1) then
      flag = NODATASETRANK
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName,&
      & 'Failed to find rank of data set.')
    endif

    ! allocate, based on above SD, dim info
    call allocate_test ( edge,   rank, 'edge',   ModuleName )
    call allocate_test ( start,  rank, 'start',  ModuleName )
    call allocate_test ( stride, rank, 'stride', ModuleName )

    ! Set "slab" dimensions
    edge = dim_sizes(1:rank)
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

    if ( present ( firstMAF ) ) then
      if ( (firstMAF >= numMAFs) .or. (firstMAF < 0) ) then
        flag = FIRSTMAFNOTFOUND
        if ( MyNeverFail ) return
        call MLSMEssage ( MLSMSG_Error, ModuleName, &
        & input_err // 'firstMAF' )
      endif
      l1bData%firstMAF = firstMAF
    else
      l1bData%firstMAF = 0
    end if

    if ( present (lastMAF) ) then
      if ( lastMAF < l1bData%firstMAF ) then
        flag = LASTMAFNOTFOUND
        if ( MyNeverFail ) return
        call MLSMEssage ( MLSMSG_Error, ModuleName, &
        & input_err // 'last' )
      endif
      if ( lastMAF >= numMAFs ) then
        l1bData%noMAFs = numMAFs - l1bData%firstMAF
      else
        l1bData%noMAFs = lastMAF - l1bData%firstMAF + 1
      end if
    else
      l1bData%noMAFs = numMAFs - l1bData%firstMAF
    end if

    noMAFs = l1bData%noMAFs
    edge(rank) = l1bData%noMAFs
    start(rank) = l1bData%firstMAF

    ! allocate, read counterMAF
    call Allocate_test ( l1bData%counterMaf, l1bData%noMAFs, &
      & 'counterMAF', ModuleName )
    status = sfrdata_f90(sds1_id,  (/ l1bData%firstMAF /) , (/1/), &
      & (/l1bData%noMAFs/), l1bData%counterMAF )
    if ( status == -1 ) then
      flag = CANTREADCOUNTERMAF
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_L1BRead // 'counterMAF.' )
    endif

    ! allocate, read according to field type; nullify unused pointers
    select case ( data_type )

    case ( DFNT_CHAR8 ) ! ----------------------- character
      allocate ( l1bData%charField(l1bData%noAuxInds,l1bData%maxMIFs, &
        & l1bData%noMAFs), STAT=alloc_err )
      if ( alloc_err /= 0 ) then
        flag = CANTALLOCATECHARS
        if ( MyNeverFail ) return
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_allocate // ' charField pointer.' )
      endif
      nullify (l1bData%intField, l1bData%dpField)
      status = sfrcdata ( sds2_id, start, stride, edge, &
        & l1bData%charField )

    case ( DFNT_INT32 ) ! ------------------------------ integer
      call allocate_test ( l1bData%intField, &
        & l1bData%noAuxInds, l1bData%maxMIFs, l1bData%noMAFs, &
        & 'l1bData%intField', ModuleName )
      nullify ( l1bData%charField, l1bData%dpField )
      status = sfrdata_f90 ( sds2_id, start, stride, edge, &
          l1bData%intField )

    case ( DFNT_FLOAT64 ) ! ------------------------------ real (r8)
      call allocate_test ( l1bData%dpField, &
        & l1bData%noAuxInds, l1bData%maxMIFs, l1bData%noMAFs, &
        & 'l1bData%dpField', ModuleName )
      nullify ( l1bData%charField, l1bData%intField )
      status = sfrdata_f90 ( sds2_id, start, stride, edge, &
          l1bData%dpField )

    case ( DFNT_FLOAT32 ) ! ------------------------------ real (r4)
      call allocate_test ( l1bData%dpField, &
        & l1bData%noAuxInds, l1bData%maxMIFs, l1bData%noMAFs, &
        & 'l1bData%dpField', ModuleName )
      nullify ( l1bData%charField, l1bData%intField )

      call allocate_test ( tmpr4Field, &
        & l1bData%noAuxInds, l1bData%maxMIFs, l1bData%noMAFs, &
        & 'tmpR4Field', ModuleName )
      status = sfrdata_f90 ( sds2_id, start, stride, edge, &
          tmpR4Field )
      l1bData%dpField = tmpR4Field
      
      call deallocate_test ( tmpr4Field, 'tmpr4Field', ModuleName )

    end select ! ----------------------------------------------------

    if ( status == -1 ) then
      flag = CANTREAD3DFIELD
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_L1BRead // quantityName )
    endif

    ! Terminate access to the data sets
    status = sfendacc(sds1_id)
    if ( status == -1 ) then
      flag = CANTENDCOUNTERMAF
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Failed to terminate access to data sets.' )
      flag = -1
    end if

    status = sfendacc(sds2_id)
    if ( status == -1 ) then
      flag = CANTENDQUANTITY
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Failed to terminate access to data sets.' )
      flag = -1
    end if

    ! deallocate local variables
    call Deallocate_test ( edge,   'edge',   ModuleName )
    call Deallocate_test ( start,  'start',  ModuleName )
    call Deallocate_test ( stride, 'stride', ModuleName )

  end subroutine ReadL1BData

  !---------------------------------------------------- DeallocateL1BData ----
  subroutine DeallocateL1BData ( l1bData )
    ! This should be called when an l1bData is finished with
    type( L1BData_T ), intent(inout) :: L1bData

    ! Local variables
    integer :: DEALLOC_ERR              ! Flag

    ! Executable code
    if ( associated(l1bData%counterMAF) ) &
      call deallocate_test ( l1bData%counterMAF,&
      & 'l1bData%counterMAF', ModuleName )
    if ( associated(l1bData%charField) ) then
      deallocate ( l1bData%charField, STAT=dealloc_err )
      if ( dealloc_err /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'l1bData%charField' )
    end if
    if ( associated(l1bData%intField) ) &
      call deallocate_test ( l1bData%intField,&
      & 'l1bData%intField', ModuleName )
    if ( associated(l1bData%dpField) ) &
      call deallocate_test ( l1bData%dpField,&
      & 'l1bData%dpField', ModuleName )
  end subroutine DeallocateL1BData

  !---------------------------------------------------- DumpL1BData ----
  subroutine DumpL1BData ( l1bData, details )
    ! Disclose pertinent, perhaps damning facts about an l1brad quantity
    type( L1BData_T ), intent(inout) :: L1bData
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even counterMAF
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 1

    ! Local variables
    integer :: MYDETAILS

    ! Executable code
    myDetails = 1
    if ( present(details) ) myDetails = details
    
    if ( myDetails < -1 ) return
    call output('L1B rad quantity Name = ', advance='no')
    call output(L1bData%L1BName, advance='yes')
    call output('  First major frame read = ', advance='no')
    call output(L1bData%FirstMAF, advance='yes')
    call output('  Num of MAFs read = ', advance='no')
    call output(L1bData%NoMAFs, advance='yes')
    call output('  Max # of MIFs/MAF in SD array = ', advance='no')
    call output(L1bData%MaxMIFs, advance='yes')
    call output('  Num of auxilliary indices = ', advance='no')
    call output(L1bData%NoAuxInds, advance='yes')

    if ( myDetails < 0 ) return
    if ( associated(l1bData%counterMAF) ) then
      call dump ( l1bData%counterMAF,&
      & 'l1bData%counterMAF' )
    else
      call output('(CounterMAF array unassociated)', advance='yes')
    endif

    if ( myDetails < 1 ) return
    if ( associated(l1bData%charField) ) then
      call dump ( l1bData%CharField,&
      & 'l1bData%CharField' )
    else
      call output('(CharField array unassociated)', advance='yes')
    endif

    if ( associated(l1bData%intField) ) then
      call dump ( l1bData%intField,&
      & 'l1bData%intField' )
    else
      call output('(intField array unassociated)', advance='yes')
    endif

    if ( associated(l1bData%dpField) ) then
      call dump ( l1bData%dpField,&
      & 'l1bData%dpField' )
    else
      call output('(dpField array unassociated)', advance='yes')
    endif
  end subroutine DumpL1BData

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
  end subroutine announce_error

end module L1BData

! $Log$
! Revision 2.14  2001/10/26 23:11:10  pwagner
! Provides a single dump module interface
!
! Revision 2.13  2001/10/25 23:31:28  pwagner
! Fixed dump; readl1bData takes neverfail option
!
! Revision 2.12  2001/10/23 22:44:15  pwagner
! Added DumpL1BData
!
! Revision 2.11  2001/10/23 17:09:04  pwagner
! Added PRECISIONSUFFIX as parameter
!
! Revision 2.10  2001/10/03 22:50:03  vsnyder
! Add unfound-quantity name to an error message
!
! Revision 2.9  2001/06/01 02:06:45  livesey
! Bug fix with counterMAF
!
! Revision 2.8  2001/05/30 23:51:48  livesey
! New version, uses new HDF, also cleaner
!
! Revision 2.7  2001/05/06 20:53:47  pwagner
! Allocates l1binfo%filenames along with ids
!
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
