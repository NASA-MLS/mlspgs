
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L1BData
!===============================================================================

  use Hdf, only: DFACC_READ, SFSTART
   USE MLSCommon, only: R8, L1BInfo_T, FileNameLen
   USE MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error, &
     & MLSMSG_L1Bread, MLSMSG_Warning
   use Lexer_Core, only: PRINT_SOURCE
  use MoreTree, only: Get_Boolean, Get_Field_ID
  use Output_M, only: Output
  use String_Table, only: Get_String
  use TREE, only: NSONS, &
    &             SUB_ROSA, SUBTREE, DUMP_TREE_NODE, SOURCE_REF
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Contents:

! Definition -- L1BData_T
! Subroutines -- ReadL1BData
!                DeallocateL1BData
!                FindL1BData

! Remarks:  This is a prototype module containing parameters, a derived type
!           definition, and a subroutine for reading the L1B data.

! Parameters

   INTEGER, PARAMETER :: NAME_LEN = 64		! Max len of SDS array name

! This data type is used to store quantities from an L1B data file.

   TYPE L1BData_T
      CHARACTER (LEN=NAME_LEN) :: L1BName	! Name of field in file
      INTEGER :: firstMAF   			! First major frame read
      INTEGER :: noMAFs				! # of MAFs read
      INTEGER :: maxMIFs			! Max # of MIFs/MAF in SD array
      INTEGER :: noAuxInds			! # of auxilliary indices

      INTEGER, DIMENSION(:), POINTER :: counterMAF	! dimensioned (noMAFs)

      CHARACTER, DIMENSION(:,:,:), POINTER :: charField
      REAL(r8), DIMENSION(:,:,:), POINTER :: dpField
      INTEGER, DIMENSION(:,:,:), POINTER :: intField
			! dimensioned (noAuxInds,maxMIFs,noMAFs)
   END TYPE L1BData_T

  private ::  announce_error
  integer, private :: Error            ! Error level -- 0 = OK

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE l1bradSetup ( root, l1bInfo, F_FILE, &
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
   END SUBROUTINE l1bradSetup
!----------------------------

!------------------------------------------------------------------------------
   SUBROUTINE l1boaSetup ( root, l1bInfo, F_FILE )
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
   END SUBROUTINE l1boaSetup
!----------------------------

!------------------------------------------------------------------------------
   SUBROUTINE ReadL1BData (L1FileHandle, quantityName, l1bData, noMAFs, flag, &
                           firstMAF, lastMAF)
!------------------------------------------------------------------------------

! Brief description of subroutine
! This is a prototype for the ReadL1BData subroutine in the Open/Init task.

! Arguments

      CHARACTER (LEN=*), INTENT(IN) :: quantityName

      INTEGER, INTENT(IN) :: L1FileHandle

      INTEGER, INTENT(IN), OPTIONAL :: firstMAF, lastMAF
   
      TYPE( L1BData_T ), INTENT(OUT) :: l1bData

      INTEGER, INTENT(OUT) :: flag, noMAFs

! Parameters

      CHARACTER (LEN=*), PARAMETER :: INPUT_ERR = 'Error in input argument '

      INTEGER, PARAMETER :: MAX_VAR_DIMS = 32

! Functions

      INTEGER, EXTERNAL :: sfginfo, sfn2index, sfselect, sfrdata, sfrcdata
      INTEGER, EXTERNAL :: sfendacc

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: alloc_err, data_type, dealloc_err, i, n_attrs, numMAFs, rank
      INTEGER :: sds1_id, sds2_id, sds_index, status
      INTEGER :: dim_sizes(MAX_VAR_DIMS)
      INTEGER, ALLOCATABLE :: edge(:), start(:), stride(:)

      LOGICAL :: firstCheck, lastCheck

      REAL, ALLOCATABLE :: realField(:,:,:)

      flag = 0

! Find data sets for counterMAF & quantity by name

      sds_index = sfn2index(L1FileHandle, 'counterMAF')
      IF (sds_index == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed &
                                       &to find index of counterMAF data set.')

      sds1_id = sfselect(L1FileHandle, sds_index)
      IF (sds1_id == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                     &find identifier of counterMAF data set.')

      sds_index = sfn2index(L1FileHandle, quantityName)
      IF (sds_index == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed &
                                         &to find index of quantity data set.')

      sds2_id = sfselect(L1FileHandle, sds_index)
      IF (sds2_id == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                             &find identifier of data set matching the index.')

! Find rank (# of dimensions), dimension sizes of quantity data set

      status = sfginfo(sds2_id, quantityName, rank, dim_sizes, data_type, &
                       n_attrs)

      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                                      &find rank of data set.')

! Allocate, based on above SD, dim info

      ALLOCATE (edge(rank), start(rank), stride(rank), STAT=alloc_err)

      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  slab dimensions.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Set "slab" dimensions

      DO i = 1, rank
         edge(i) = dim_sizes(i)
      ENDDO
      edge(rank) = 1

      start = 0

      stride = 1

! Fill in "indexing" values of l1b object

      l1bData%L1BName = quantityName

      IF (rank < 2) THEN
         l1bData%maxMIFs = 1
      ELSE IF (rank > 2) THEN
         l1bData%maxMIFs = dim_sizes(2)
      ELSE
         l1bData%maxMIFs = dim_sizes(1)
      ENDIF

      IF (rank > 2) THEN
         l1bData%noAuxInds = dim_sizes(1)
      ELSE
         l1bData%noAuxInds = 1
      ENDIF

! Check input arguments, set noMAFs

      numMAFs = dim_sizes(rank)

      firstCheck = PRESENT(firstMAF)
      lastCheck = PRESENT(lastMAF)

      IF (firstCheck) THEN

          IF ( (firstMAF >= numMAFs) .OR. (firstMAF < 0) ) THEN
             msr = INPUT_ERR // 'firstMAF'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

          l1bData%firstMAF = firstMAF

      ELSE

          l1bData%firstMAF = 0

      ENDIF

      IF (lastCheck) THEN

          IF (lastMAF < l1bData%firstMAF) THEN
             msr = INPUT_ERR // 'lastMAF'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

          IF (lastMAF >= numMAFs) THEN
             l1bData%noMAFs = numMAFs - l1bData%firstMAF
          ELSE
             l1bData%noMAFs = lastMAF - l1bData%firstMAF + 1
          ENDIF

      ELSE

         l1bData%noMAFs = numMAFs - l1bData%firstMAF

      ENDIF

      noMAFs = l1bData%noMAFs

! Allocate, read counterMAF

      ALLOCATE(l1bData%counterMAF(l1bData%noMAFs), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  counterMAF pointer.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      status = sfrdata(sds1_id,  l1bData%firstMAF , stride, &
        & l1bData%noMAFs, l1bData%counterMAF )
      IF (status == -1) THEN
         msr = MLSMSG_L1BRead // 'counterMAF.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Allocate, read according to field type; nullify unused pointers

      IF (data_type == 4) THEN

! character (DFNT_CHAR8)

         ALLOCATE ( l1bData%charField(l1bData%noAuxInds,l1bData%maxMIFs,&
                   &l1bData%noMAFs), STAT=alloc_err )
         IF ( alloc_err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' charField pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         NULLIFY (l1bData%intField, l1bData%dpField)

         DO i = 1, l1bData%noMAFs
            start(rank) = l1bData%firstMAF + (i-1)
            status = sfrcdata( sds2_id, start, stride, edge, &
                               l1bData%charField(:,:,i) )
            IF (status == -1) THEN
               msr = MLSMSG_L1BRead // quantityName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
         ENDDO

      ELSE IF (data_type == 24) THEN

! integer (DFNT_INT32)

         ALLOCATE (l1bData%intField(l1bData%noAuxInds,l1bData%maxMIFs,&
                  &l1bData%noMAFs), STAT=alloc_err)
         IF ( alloc_err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' intField pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         NULLIFY (l1bData%charField, l1bData%dpField)

         DO i = 1, l1bData%noMAFs
            start(rank) = l1bData%firstMAF + (i-1)
            status = sfrdata( sds2_id, start, stride, edge, &
                              l1bData%intField(:,:,i) )
            IF (status == -1) THEN
               msr = MLSMSG_L1BRead // quantityName
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF
         ENDDO

      ELSE

! float (REAL or DOUBLE PRECISION)

         ALLOCATE (l1bData%dpField(l1bData%noAuxInds,l1bData%maxMIFs,&
                  &l1bData%noMAFs), STAT=alloc_err)
         IF ( alloc_err /= 0 ) THEN
            msr = MLSMSG_Allocate // ' dpField pointer.'
            CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
         ENDIF
         NULLIFY (l1bData%charField, l1bData%intField)

! If real (DFNT_FLOAT32), allocate real variable and then move to dp pointer

         IF (data_type == 5) THEN

            ALLOCATE (realField(l1bData%noAuxInds,l1bData%maxMIFs,&
                     &l1bData%noMAFs), STAT=alloc_err)
            IF ( alloc_err /= 0 ) THEN
               msr = MLSMSG_Allocate // ' realField pointer.'
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

            DO i = 1, l1bData%noMAFs
               start(rank) = l1bData%firstMAF + (i-1)
               status = sfrdata( sds2_id, start, stride, edge, &
                                realField(:,:,i) )
               IF (status == -1) THEN
                  msr = MLSMSG_L1BRead // quantityName
                  CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
               ENDIF
            ENDDO

            l1bData%dpField = DBLE(realField)

            DEALLOCATE(realField, STAT=dealloc_err)
            IF ( dealloc_err /= 0 ) THEN
               CALL MLSMessage(MLSMSG_Warning, ModuleName, 'Failed &
                                     &deallocation of local real variable.')
               flag = -1
            ENDIF


         ELSE

! double precision (DFNT_FLOAT64)


            DO i = 1, l1bData%noMAFs
               start(rank) = l1bData%firstMAF + (i-1)
               status = sfrdata(sds2_id, start, stride, edge, &
                                l1bData%dpField(:,:,i))
               IF (status == -1) THEN
                  msr = MLSMSG_L1BRead // quantityName
                  CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
               ENDIF
            ENDDO

         ENDIF

      ENDIF

! Terminate access to the data sets

      status = sfendacc(sds1_id)
      status = sfendacc(sds2_id)
      IF (status == -1) THEN
         CALL MLSMessage(MLSMSG_Warning, ModuleName, 'Failed to terminate &
                                                     &access to data sets.')
         flag = -1
      ENDIF

! Deallocate local variables

      DEALLOCATE (edge, start, stride, STAT=dealloc_err )

      IF ( dealloc_err /= 0 ) THEN
         CALL MLSMessage(MLSMSG_Warning, ModuleName, 'Failed deallocation of &
                                                     &local variables.')
         flag = -1
      ENDIF

!----------------------------
   END SUBROUTINE ReadL1BData
!----------------------------

!----------------------------------------------
   SUBROUTINE DeallocateL1BData (l1bData, flag)
!----------------------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal field pointers of the L1BData_T
! derived type, after the calling program has finished with the data.

! Arguments

      TYPE( L1BData_T ), intent(INOUT) :: l1bData

      INTEGER, INTENT(OUT) :: flag

! Parameters

! Functions

! Variables

      INTEGER :: dealloc_err

      flag = 0

      IF ( ASSOCIATED(l1bData%counterMAF) ) DEALLOCATE (l1bData%counterMAF, &
                                                        STAT=dealloc_err)
      IF ( ASSOCIATED(l1bData%charField) ) THEN
         DEALLOCATE (l1bData%charField, STAT=dealloc_err)
      ELSE IF ( ASSOCIATED(l1bData%intField) ) THEN
         DEALLOCATE (l1bData%intField, STAT=dealloc_err)
      ELSE
         DEALLOCATE (l1bData%dpField, STAT=dealloc_err)
      ENDIF

      IF ( dealloc_err /= 0 ) THEN
         CALL MLSMessage(MLSMSG_Warning, ModuleName, 'Failed deallocation of &
                                                     &L1BData pointers.')
         flag = -1
      ENDIF

!----------------------------------
   END SUBROUTINE DeallocateL1BData
!----------------------------------

   ! ------------------------------------------- FindL1BData ----
   integer function FindL1BData( files, fieldName )
     integer, dimension(:), intent(in) :: files ! File handles
     character (len=*), intent(in) :: fieldName ! Name of field
     
     ! Externals
     integer, external :: SFN2INDEX

     ! Local variables
     integer :: i

     ! Executable code
     FindL1BData=0
     do i=1,size(files)
       if ( sfn2index(files(i),fieldName) /= -1) then
         FindL1BData=files(i)
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
END MODULE L1BData
!=================

! $Log$
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
