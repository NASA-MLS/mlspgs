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
  use L2AUXData, only: L2AUXDATA_T, WriteL2AUXData
  use L2GPData, only: L2GPData_T, WriteL2GPData, L2GPNameLen
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: I4
  use MLSFiles, only: GetPCFromRef
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSPCF2, only: MLSPCF_L2DGM_END, MLSPCF_L2DGM_START, MLSPCF_L2GP_END, &
    & MLSPCF_L2GP_START, &
	 & mlspcf_mcf_l2gp_start, mlspcf_mcf_l2dgm_start, mlspcf_mcf_l2log_start, &
	 & mlspcf_mcf_l2dgg_start
  use MoreTree, only: Get_Spec_ID
  use OUTPUT_M, only: OUTPUT
  use SDPToolkit, only: Pgs_pc_getReference, PGS_S_SUCCESS, Pgs_smf_getMsg
  use STRING_TABLE, only: GET_STRING
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TOGGLES, only: GEN, TOGGLE
  use TREE, only: DECORATION, DUMP_TREE_NODE, NODE_ID, NSONS, SOURCE_REF, &
    & SUBTREE, SUB_ROSA
  use TREE_TYPES, only: N_NAMED
  use WriteMetadata, only: PCFData_T, populate_metadata_std, &
  & populate_metadata_oth, WriteMetaLog, get_l2gp_mcf

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

contains ! =====     Public Procedures     =============================

  ! -----------------------------------------------  Output_Close  -----
  subroutine Output_Close ( root, l2gpDatabase, l2auxDatabase, l2pcf, anText )

	! Hard-wired assumptions:
	
	! ----------------------- metadata ------------------------
	
	!   for the l2aux the mcf is mlspcf_mcf_l2dgm_start
	!   for the log file the mcf is mlspcf_mcf_l2log_start
	!   for the dgg file the mcf is mlspcf_mcf_l2dgg_start

	! The correspondence between MCF and l2gp files is determined by
	! the value of        MCFFORL2GPOPTION
	! Either
	!                          (1)
	! The PCF numbers for the mcf corresponding to each
	! of the l2gp files begin with mlspcf_mcf_l2gp_start
	! and increase 1 by 1 with each succeeding species.
	! Then, after the last single-species l2gp, the very next pcf number
	! is for the one called 'other' ML2OTH.001.MCF
	! This hateful inflexibility leads to possibility
	
	!                          (2)
	! Each l2gp file name, stripped of their paths, fits the pattern
	!  *_l2gp_species_*
	! and the corresponding MCF files fit the pattern
	!  *SPECIES.*
	! where species and SPECIES are case-nsensitive "species" name
	! i.e., BrO, ClO, etc.
	! Warning:
	! You therefore must use exactly the same abbreviation for the l2gp and the
	! corresponding MCF: if the MCF is ML2T.001.MCF, don't use "temp"
	! in the l2gp name
	! This inflexibility replaces the different kind in option (1)

  ! Arguments
    integer, intent(in) :: ROOT   ! Of the output section's AST
    type (L2GPData_T), dimension(:), pointer :: l2gpDatabase ! L2GP products
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase ! L2AUX products
	 type(PCFData_T) :: l2pcf
      CHARACTER (LEN=1), POINTER :: anText(:)


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
    integer :: l2gp_mcf, l2aux_mcf, l2dgg_mcf  ! mcf numbers for writing metadata
    character (len=132) :: l2auxPhysicalFilename
    integer :: l2gpFileHandle, l2gp_Version
    character (len=132) :: l2gpPhysicalFilename
    integer, parameter:: MAXQUANTITIESPERFILE=64        
    integer, parameter :: MCFFORL2GPOPTION=1		! Either 1 or 2
    character (len=32) :: mnemonic
    character (len=256) :: msg
    integer :: NAME                     ! string index of label on output
    integer:: numquantitiesperfile        
    integer :: OUTPUT_TYPE              ! L_L2AUX or L_L2GP
    character(len=L2GPNameLen), dimension(MAXQUANTITIESPERFILE) :: QuantityNames  ! From "quantities" field
    integer :: returnStatus
    integer(i4) :: SDFID                ! File handle
    integer :: SON                      ! Of Root -- spec_args or named node
    integer :: SPEC_NO                  ! Index of son of Root
    integer(i4) :: START(3), STRIDE(3), DIM_SIZES(3)
    integer :: SWFID
    REAL :: T1, T2     ! for timing
    logical :: TIMING

    ! Executable code
    timing = .false.

    if ( toggle(gen) ) call trace_begin ( "Output_Close", root)

    error = 0
    got = .false.
    l2gp_Version = 1
    l2aux_Version = 1

    l2gp_mcf = mlspcf_mcf_l2gp_start
    l2aux_mcf = mlspcf_mcf_l2dgm_start
    l2dgg_mcf = mlspcf_mcf_l2dgg_start

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

      select case( get_spec_id(key) )
      case ( s_output )
        do field_no = 2, nsons(key)       ! Skip the command name
          gson = subtree(field_no, key)   ! An assign node
          field_index = decoration(subtree(1,gson))
          got(field_index) = .true.
          select case ( field_index )   ! Field name
          case ( f_file )
            call get_string ( sub_rosa(subtree(2,gson)), file_base )
            file_base=file_base(2:LEN_TRIM(file_base)-1) ! Parser includes quotes
          case ( f_type )
            output_type = decoration(subtree(2,gson))
          case default                  ! Everything else processed later
          end select
        end do

        select case ( output_type )
        case ( l_l2gp )

          ! Get the l2gp file name from the PCF

 !         l2gpFileHandle = mlspcf_l2gp_start
 !         found = .FALSE.
 !         do while ((.NOT. found) .AND. (l2gpFileHandle <=  mlspcf_l2gp_end))
 !           l2gp_Version=1
 !           returnStatus = Pgs_pc_getReference(l2gpFileHandle, l2gp_Version, &
 !             & l2gpPhysicalFilename)
 !           
 !           if ( returnStatus == PGS_S_SUCCESS ) then
 !             if ( INDEX(l2gpPhysicalFilename, TRIM(file_base)) /= 0 ) then
 !               found = .true.
 
  				l2gpFileHandle = GetPCFromRef(file_base, mlspcf_l2gp_start, &
				& mlspcf_l2gp_end, &
  & .TRUE., returnStatus, l2gp_Version)
 
            if ( returnStatus == 0 ) then
                ! Open the HDF-EOS file and write swath data
                
                swfid = swopen(file_base, DFACC_CREATE)

                ! Loop over the segments of the l2cf line
                
					numquantitiesperfile = 0
                do field_no = 2, nsons(key) ! Skip "output" name
                  gson = subtree(field_no,key)
                  select case ( decoration(subtree(1,gson)) )
                  case ( f_quantities )
                    do in_field_no = 2, nsons(gson)
                      db_index = -decoration(decoration(subtree(in_field_no ,gson)))
                      CALL WriteL2GPData(l2gpDatabase(db_index),swfid)
							 numquantitiesperfile=numquantitiesperfile+1
							 QuantityNames(numquantitiesperfile) = l2gpDatabase(db_index)%name
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

					! Write the metadata file
					if(MCFFORL2GPOPTION == 2) then
						l2gp_mcf = get_l2gp_mcf(swfid)
					endif
					if(l2gp_mcf <= 0) then
						call announce_error(son, &
						& 'No mcf numbers correspond to this l2gp file', l2gp_mcf)
					elseif(numquantitiesperfile <= 0) then
						call announce_error(son, &
						& 'No quantities written for this l2gp file')
					elseif(QuantityNames(numquantitiesperfile) &
					& == QuantityNames(1) ) then

					! Typical homogeneous l2gp file: e.g., BrO is ML2BRO.001.MCF
						call populate_metadata_std &
						& (swfid, l2gp_mcf, l2pcf, QuantityNames(1), anText)
						l2gp_mcf = l2gp_mcf + 1

					else

					! Type l2gp file 'other'
						call populate_metadata_oth &
						& (swfid, l2gp_mcf, l2pcf, &
						& numquantitiesperfile, QuantityNames, anText)

					endif
!              else                ! Found the right file
!                l2gpFileHandle = l2gpFileHandle + 1
!              end if

            else
!              call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
!              call MLSMessage ( MLSMSG_Error, ModuleName, &
!                &  "Error finding l2gp file:  "//mnemonic//" "//msg )
              call announce_error ( ROOT, &
                &  "Error finding l2gp file:  "//file_base, returnStatus)

            end if
            
!          end do !(.not. found) .and. (l2gpFileHandle <=  mlspcf_l2gp_end)
!          IF (.NOT. found) CALL MLSMessage(MLSMSG_Error,ModuleName,&
!            'Unable to find filename containing '//TRIM(file_base))
          
        case ( l_l2aux )
          
          ! Get the l2aux file name from the PCF
          
!          l2auxFileHandle = mlspcf_l2dgm_start
!          found = .FALSE.
!          do while ((.NOT. found) .AND. (l2auxFileHandle <=  mlspcf_l2dgm_end))
!            l2aux_Version=1
!            returnStatus = Pgs_pc_getReference(l2auxFileHandle, l2aux_Version, &
!              & l2auxPhysicalFilename)
            
!            if ( returnStatus == PGS_S_SUCCESS ) then
!              if ( INDEX(l2auxPhysicalFilename, TRIM(file_base)) /= 0 )then
!                found = .TRUE.
                
  				l2auxFileHandle = GetPCFromRef(file_base, mlspcf_l2dgm_start, &
				& mlspcf_l2dgm_end, &
  & .TRUE., returnStatus, l2aux_Version)
 
            if ( returnStatus == 0 ) then
				
                ! Create the HDF file and initialize the SD interface
                sdfId = sfstart(l2auxPhysicalFilename, DFACC_CREATE)
                
					numquantitiesperfile = 0
                do field_no = 2, nsons(key) ! Skip "output" name
                  gson = subtree(field_no,key)
                  select case ( decoration(subtree(1,gson)) )
                  case ( f_quantities )
                    do in_field_no = 2, nsons(gson)
                      db_index = -decoration(decoration(subtree(in_field_no ,gson)))
                      CALL WriteL2AUXData(l2auxDatabase(db_index),sdfid)
							 numquantitiesperfile=numquantitiesperfile+1
			            call get_string ( sub_rosa(subtree(2,gson)), QuantityNames(numquantitiesperfile) )
                    end do ! in_field_no = 2, nsons(gson)
                  case ( f_overlaps )
                    ! ??? More work needed here
                  end select
                end do ! field_no = 2, nsons(key)
                
                ! Now close the file
                returnStatus = sfend(sdfid)
                if ( returnStatus /= PGS_S_SUCCESS ) then
                  call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
                  call MLSMessage ( MLSMSG_Error, ModuleName, &
                    &  "Error closing l2aux file:  "//mnemonic//" "//msg )
                end if

				! Write the metadata file
					if(numquantitiesperfile <= 0) then
						call announce_error(son, &
						& 'No quantities written for this l2aux file')
					else

					! We will need to think harder about this; until then reuse
						call populate_metadata_oth &
						& (swfid, l2aux_mcf, l2pcf, &
						& numquantitiesperfile, QuantityNames, anText)

					endif
!              else
!                l2auxFileHandle =  l2auxFileHandle + 1
!              end if

            else
!              call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
!              call MLSMessage ( MLSMSG_Error, ModuleName, &
!                &  "Error finding l2aux file:  "//mnemonic//" "//msg )
              call announce_error ( ROOT, &
                &  "Error finding l2aux file:  "//file_base, returnStatus)
            end if
            
!          end do ! (.not. found) .and. (l2auxFileHandle <=  mlspcf_l2aux_end)

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

    if ( toggle(gen) ) call trace_end ( "Output_Close")

! Write the log file metadata

      CALL WriteMetaLog(l2pcf)


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
  subroutine ANNOUNCE_ERROR ( WHERE, full_message, CODE )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
	character(LEN=*), intent(in)    :: full_message
    integer, intent(in), optional :: CODE    ! Code for error message

    error = max(error,1)
    call output ( '***** At ' )
!    call print_source ( source_ref(where) )
	if(where > 0) then
	    call print_source ( source_ref(where) )
		else
    call output ( '(no lcf node available)' )
		endif
    call output ( ' OutputAndClose complained: ' )


		CALL output("Caused the following error:", advance='yes', &
		& from_where=ModuleName)
		CALL output(trim(full_message), advance='yes', &
		& from_where=ModuleName)
		if(present(code)) then
			select case ( code )
			end select
		endif
    end subroutine ANNOUNCE_ERROR

end module OutputAndClose

! $Log$
! Revision 2.15  2001/04/03 23:51:28  pwagner
! Many changes; some may be right
!
! Revision 2.14  2001/04/02 23:43:46  pwagner
! Now makes metadata calls; it compiles, but does it bomb?
!
! Revision 2.13  2001/03/28 00:23:20  pwagner
! Made tiny changes to use announce_error
!
! Revision 2.12  2001/03/20 18:35:02  pwagner
! Using GetPCFromRef to get file handles
!
! Revision 2.11  2001/03/15 21:18:57  vsnyder
! Use Get_Spec_ID instead of decoration(subtree...
!
! Revision 2.10  2001/03/06 22:40:24  livesey
! Working l2aux
!
! Revision 2.9  2001/02/23 18:15:48  livesey
! Added trace calls.
!
! Revision 2.8  2001/02/13 22:59:36  pwagner
! l2 modules can only use MLSPCF2
!
! Revision 2.7  2001/02/09 19:30:16  vsnyder
! Move checking for required and duplicate fields to init_tables_module
!
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

