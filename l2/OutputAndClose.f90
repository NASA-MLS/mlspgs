! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=======================================================================================

module OutputAndClose ! outputs all data from the Join module to the
                      ! appropriate L2 Files

!=======================================================================================

  use Hdf, only: DFACC_CREATE, DFNT_FLOAT64, SFCREATE, SFDIMID, SFEND, &
    & SFENDACC, SFSTART, SFWDATA
  use HDFEOS, only: SWCLOSE, SWOPEN
  use INIT_TABLES_MODULE, only: F_FILE, F_OVERLAPS, F_QUANTITIES, F_TYPE, &
    & FIELD_FIRST, FIELD_LAST, L_L2AUX, L_L2GP, S_OUTPUT, S_TIME
  use L2AUXData, only: L2AUXDATA_T, L2AUXDIMNAMES
  use L2GPData, only: L2GPData_T, WriteL2GPData
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: I4
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSPCF, only: MLSPCF_L2AUX_END, MLSPCF_L2AUX_START, MLSPCF_L2GP_END, &
    & MLSPCF_L2GP_START
  use OUTPUT_M, only: OUTPUT
  use SDPToolkit, only: Pgs_pc_getReference, PGS_S_SUCCESS, Pgs_smf_getMsg
  use STRING_TABLE, only: GET_STRING
  use TREE, only: DECORATION, DUMP_TREE_NODE, NODE_ID, NSONS, SOURCE_REF, &
    & SUBTREE, SUB_ROSA
  use TREE_TYPES, only: N_NAMED

  implicit none
  private
  public :: Output_Close

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130), private :: id = & 
    & "$Id$"
  character(len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  !-----------------------------------------------------------------------------

  ! -----     Private declarations     ---------------------------------

  ! For Announce_Error
    integer :: ERROR
    integer, parameter :: DUPLICATE_FIELD = 1

contains ! =====     Public Procedures     =============================

  ! -----------------------------------------------  Output_Close  -----
  subroutine Output_Close ( root, l2gpDatabase, l2auxDatabase )

  ! Arguments
    integer, intent(in) :: ROOT   ! Of the output section's AST
    type (L2GPData_T), dimension(:), pointer :: l2gpDatabase ! L2GP products
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase ! L2AUX products

  ! - - - External Procedures

    interface
      integer function SFSDMNAME ( DIM_ID, DIM_NAME ) ! An HDF function
        integer, intent(in) :: DIM_ID
        character(len=*), intent(in) :: DIM_NAME
      end function SFSDMNAME
    end interface

  ! - - - Local declarations - - -

    integer :: db_index
    integer :: FIELD_INDEX              ! F_... field code
    integer :: FIELD_NO                 ! Index of assign vertex sons of Key
    character (len=132) :: FILE_BASE    ! From the FILE= field
    integer :: FLAG
    logical :: FOUND
    logical :: GOT(field_first:field_last)
    integer :: GSON                     ! Son of Son -- an assign node
    integer :: I
    integer :: IN_FIELD                 ! A son of an assign vertex
    integer :: IN_FIELD_NO              ! Index of sons of assign vertex
    integer :: KEY                      ! Index of spec_args node
    integer :: l2auxFileHandle, l2aux_Version
    character (len=132) :: l2auxPhysicalFilename
    integer :: l2gpFileHandle, l2gp_Version
    character (len=132) :: l2gpPhysicalFilename
    character (len=32) :: mnemonic
    character (len=256) :: msg
    integer :: NAME                     ! string index of label on output
    integer :: OUTPUT_TYPE              ! L_L2AUX or L_L2GP
    character(len=132) :: QuantityName  ! From "quantities" field
    integer :: returnStatus
    integer(i4) :: SD_ID, SDS_ID, DIM_ID
    integer :: SON                      ! Of Root -- spec_args or named node
    integer :: SPEC_NO                  ! Index of son of Root
    integer(i4) :: START(3), STRIDE(3), DIM_SIZES(3)
    integer :: SWFID
    REAL :: T1, T2     ! for timing
    logical :: TIMING

    ! Executable code
    timing = .false.

    error = 0
    got = .false.
    l2gp_Version = 1
    l2aux_Version = 1

    ! Loop over the lines in the l2cf

    do spec_no = 2, nsons(root)-1 ! Skip name at begin and end of section
      son = subtree(spec_no,root)
      if ( node_id(son) == n_named ) then ! Is spec labeled?
        key = subtree(2,son)
        name = sub_rosa(subtree(1,son))
      else ! Son is n_spec_args
        key = son
        name = 0
      end if

      select case( decoration(subtree(1,decoration(subtree(1,key)))) )
      case ( s_output )
        do field_no = 2, nsons(key)       ! Skip the command name
          gson = subtree(field_no, key)   ! An assign node
          field_index = decoration(subtree(1,gson))
          if ( got(field_index) ) then
            call announce_error ( gson, duplicate_field )
          else
            got(field_index) = .true.
            select case ( field_index )   ! Field name
            case ( f_file )
              call get_string ( sub_rosa(subtree(2,gson)), file_base )
              file_base=file_base(2:LEN_TRIM(file_base)-1) ! Parser includes quotes
            case ( f_type )
              output_type = decoration(subtree(2,gson))
            case default                  ! Everything else processed later
            end select
          end if
        end do

        select case ( output_type )
        case ( l_l2gp )

          ! Get the l2gp file name from the PCF

          l2gpFileHandle = mlspcf_l2gp_start
          found = .FALSE.
          do while ((.NOT. found) .AND. (l2gpFileHandle <=  mlspcf_l2gp_end))
            l2gp_Version=1
            returnStatus = Pgs_pc_getReference(l2gpFileHandle, l2gp_Version, &
              & l2gpPhysicalFilename)
            
            if ( returnStatus == PGS_S_SUCCESS ) then
               if ( INDEX(l2gpPhysicalFilename, TRIM(file_base)) /= 0 ) then
                  found = .true.
                  ! Open the HDF-EOS file and write swath data

                  swfid = swopen(l2gpPhysicalFilename, DFACC_CREATE)

                  ! Loop over the segments of the l2cf line

                  do field_no = 2, nsons(key) ! Skip "output" name
                     gson = subtree(field_no,key)
                     select case ( decoration(subtree(1,gson)) )
                     case ( f_quantities )
                        do in_field_no = 2, nsons(gson)
                           db_index = -decoration(decoration(subtree(in_field_no ,gson)))
                           CALL WriteL2GPData(l2gpDatabase(db_index),swfid,swathName='Fred')
                        end do ! in_field_no = 2, nsons(gson)
                     case ( f_overlaps )
                        ! ??? More work needed here
                     end select
                  end do ! field_no = 2, nsons(key)
                  returnStatus = swclose(swfid)
                  if (returnStatus /= PGS_S_SUCCESS) then
                     call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
                     call MLSMessage ( MLSMSG_Error, ModuleName, &
                          &  "Error closing  l2gp file:  "//mnemonic//" "//msg )
                  end if

               else                ! Found the right file
                  l2gpFileHandle = l2gpFileHandle + 1
               end if
            else
               call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
               call MLSMessage ( MLSMSG_Error, ModuleName, &
                    &  "Error finding  l2gp file:  "//mnemonic//" "//msg )
            end if

        end do !(.not. found) .and. (l2gpFileHandle <=  mlspcf_l2gp_end)
        IF (.NOT. found) CALL MLSMessage(MLSMSG_Error,ModuleName,&
             'Unable to find filename containing '//TRIM(file_base))

        case ( l_l2aux )

          ! Get the l2aux file name from the PCF

          l2auxFileHandle = mlspcf_l2aux_start
          found = .FALSE.
          do while ((.NOT. found) .AND. (l2auxFileHandle <=  mlspcf_l2aux_end))
             l2aux_Version=1
            returnStatus = Pgs_pc_getReference(l2auxFileHandle, l2aux_Version, &
              & l2auxPhysicalFilename)

            if ( returnStatus == PGS_S_SUCCESS ) then
              if ( INDEX(file_base, l2auxPhysicalFilename) /= 0 )then
                found = .TRUE.

                ! Create the HDF file and initialize the SD interface

                sd_id = sfstart(l2auxPhysicalFilename, DFACC_CREATE)

                ! Loop over the segments of the l2aux line

                do field_no = 2, nsons(key) ! Skip "output" name
                  gson = subtree(field_no,key) ! An assign vertex
                  select case ( decoration(subtree(1,gson)) )
                  case ( f_quantities )
                    do in_field_no = 2, nsons(gson) ! Skip "quantities" name
                      in_field = subtree(in_field_no,gson)
                      db_index = -decoration(decoration(in_field_no))

                      ! Set up dimensions
                      dim_sizes(1:3) = &
                        & l2auxDatabase(db_index)%dimensions(1:3)%noValues
                      if ( l2auxDatabase(db_index)%noDimensionsUsed == 2 ) &
                        & dim_sizes(3)=1
                      ! Create data set
                      call get_string ( sub_rosa(in_field), quantityName )
                      sds_id = sfcreate(sd_id, quantityName, DFNT_FLOAT64, &
                        & l2auxDatabase(db_index)%noDimensionsUsed, dim_sizes)

                      ! Give names to the dimensions
                      do i = 1, l2auxDatabase(db_index)%noDimensionsUsed
                        dim_id = sfdimid(sds_id, i)
                        returnStatus = sfsdmname(dim_id, &
                          & L2AUXDimNames(l2auxDatabase(db_index)% &
                          & dimensions(i)%dimensionFamily))
                        if ( returnStatus /= PGS_S_SUCCESS ) then
                          call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
                          call MLSMessage ( MLSMSG_Error, ModuleName, &
                            & "Error setting dimension name to SDS l2aux file:  "// &
                            & L2AUXDimNames(l2auxDatabase(db_index)%dimensions(i)% &
                            & dimensionFamily) //mnemonic//" "//msg )
                        end if
                      end do

                      ! Write out the SDS data
                      start=0
                      stride=1

                      returnStatus = sfwdata(sds_id, start, stride, dim_sizes, &
                        & l2auxDatabase(db_index)%values )
                      if ( returnStatus /= PGS_S_SUCCESS ) then
                        call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
                        call MLSMessage ( MLSMSG_Error, ModuleName, &
                          & "Error writing SDS data to  l2aux file:  " // &
                          & mnemonic // " " // msg )
                      end if

                      ! Terminate access to the data set

                      returnStatus = sfendacc(sds_id)
                    end do ! in_field_no = 2, nsons(gson)
                  case ( f_overlaps )
                  ! ??? More work needed here
                  end select
                end do ! field_no = 2, nsons(key)

                returnStatus = sfend(sd_id)
                if ( returnStatus /= PGS_S_SUCCESS ) then
                  call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
                  call MLSMessage ( MLSMSG_Error, ModuleName, &
                    &  "Error closing l2aux file:  "//mnemonic//" "//msg )
                end if

              else
                l2auxFileHandle =  l2auxFileHandle + 1
              end if
            else

              call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
              call MLSMessage ( MLSMSG_Error, ModuleName, &
                &  "Error finding  l2aux file:  "//mnemonic//" "//msg )
            end if

          end do ! (.not. found) .and. (l2auxFileHandle <=  mlspcf_l2aux_end)
        end select
      case ( s_time )
        if ( timing ) then
          call sayTime
        else
          call cpu_time ( t1 )
          timing = .true.
        end if
      end select
    end do  ! spec_no
    if ( timing ) call sayTime

  contains
    subroutine SayTime
      call cpu_time ( t2 )
      call output ( "Timing for Output_Close =" )
      call output ( DBLE(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

  end subroutine Output_Close

! =====     Private Procedures     =====================================

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref(where) )
    call output ( ': ' )
    select case ( code )
    case ( duplicate_field )
      call output ( "The " )
      call dump_tree_node ( where, 0 )
      call output ( " field appears more than once.", advance='yes' )
    end select
    end subroutine ANNOUNCE_ERROR

end module OutputAndClose

! $Log$
! Revision 2.6  2001/02/09 00:38:22  livesey
! Various updates
!
! Revision 2.5  2001/01/03 18:15:13  pwagner
! Changed types of t1, t2 to real
!
! Revision 2.4  2000/11/16 02:25:13  vsnyder
! Implement timing.
!
! Revision 2.3  2000/10/05 16:43:00  pwagner
! Now compiles with new L2GPData module
!
! Revision 2.2  2000/09/11 19:43:47  ahanzel
! Removed old log entries.
!
! Revision 2.1  2000/09/08 22:55:56  vsnyder
! Revised to use the tree output by the parser
!

