! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE MLSL2Options              !  Options and Settings for the MLSL2 program
!=============================================================================

  use INTRINSIC, only: L_HOURS, L_MINUTES, L_SECONDS
  use MLSFiles, only: WILDCARDHDFVERSION, HDFVERSION_4, HDFVERSION_5
  use MLSMessageModule, only: MLSMSG_Error
  use MLSPCF2, only: MLSPCF_L1B_RAD_END, MLSPCF_L1B_RAD_START

  implicit none
  public
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module simply contains initial or permanent settings. Values
  ! are chosen according to what is most suitable for the environment.
  ! For example
  ! certain settings may be appropriate during development but not
  ! for production use, i.e. sips. Therefore, this is a convenient place
  ! to hold everything that needs to be changed before delivery.
  
  ! See also MLSL2PCF, L2ParInfo.parallel, lib/toggles.switches

  ! --------------------------------------------------------------------------
  ! The following should be adjusted before delivery to sips

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  ! Set the following to TRUE before delivering level 2 to sips
  logical, parameter :: SIPS_VERSION =  .false. 

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  ! Update these lines before delivery to sips     
  ! id to print out in response to "--version" command-line option       
  character(LEN=*), dimension(2), parameter :: CURRENT_VERSION_ID = (/ &    
    & 'v2.xx swdev team              ', &       
    & 'Copyright statement omitted   '/)
     
  ! Set the following to 1 before delivering to sips;                       
  ! when set to 0, it allows program to run w/o creating metadata           
  integer            ::                         PENALTY_FOR_NO_METADATA = 0

  ! Set the following to -2 before delivering to sips;                      
  ! (its possible values and their effects on normal output:                
  ! -1          sent to stdout (via print *, '...')                         
  ! -2          sent to Log file (via MLSMessage)                           
  ! < -2        both stdout and Log file                                    
  ! > -1        Fortran 'unit=OUTPUT_PRINT_UNIT')                           
  integer            :: OUTPUT_PRINT_UNIT = -2                              

  ! Set the following to MLSMSG_Error before delivering to sips;
  ! when set higher, it allows program keep going despite errors
  ! when set lower, the program would quit even on warnings
  integer, parameter :: QUIT_ERROR_THRESHOLD = MLSMSG_Error

  ! Set the following to 2 before delivering to sips;
  ! If 0, you won't be able to distinguish normal termination
  ! from some abnormal ones (e.g. in parser) (a bad thing)
  ! if 2, status will be 2 only if run complete                             
  ! and without error (a good thing)
  integer, parameter :: NORMAL_EXIT_STATUS = 2          

  ! ---------------------------------------------------------------
  ! None of the following need to be changed before delivery to sips
  
  ! Assume hdf files w/o explicit hdfVersion field are this                 
  ! 4 corresponds to hdf4, 5 to hdf5 in L2GP, L2AUX, etc.                   
  integer            :: DEFAULT_HDFVERSION_WRITE = HDFVERSION_5
  ! Set to WILDCARDHDFVERSION if you wish to autodetect such files          
  ! on input                                                                
  integer            :: DEFAULT_HDFVERSION_READ = WILDCARDHDFVERSION
  integer            :: LEVEL1_HDFVERSION = WILDCARDHDFVERSION

  ! What units to use in summarizing timings at end of run
  integer            :: SECTIONTIMINGUNITS = L_SECONDS
  logical            :: patch = .false.       ! Set if run must not create file,
  ! Whether to restart printing identical warnings at each new phase
  logical            :: RESTARTWARNINGS = .true.
  ! Whether to skip doing the direct writes--quicker when snooping       swath  
  logical            :: SKIPDIRECTWRITES = .false.         
  logical            :: SKIPDIRECTWRITESORIGINAL = .false.         
  ! Whether to skip doing the retrieval--a pre-flight checkout of paths, etc.
  logical            :: SKIPRETRIEVAL = .false.
  logical            :: SKIPRETRIEVALORIGINAL = .false. ! May skip for some phases
  ! In case special dumps are to go to a special dumpfile
  character(len=255) :: SPECIALDUMPFILE = ' '
  ! Whether to stop after the chunk division section
  logical            :: STOPAFTERCHUNKDIVIDE = .false.         
  ! Whether to stop after doing the global settings section
  logical            :: STOPAFTERGLOBAL = .false.         
  ! Whether to exit with status 1 no matter what
  logical            :: STOPWITHERROR = .false.         
  ! Whether to do only a pre-flight checkout of paths
  logical            :: CHECKPATHS = .false.         
  ! Whether to catenate split autoDirectWrites
  logical            :: CATENATESPLITS = .false.         

  logical            :: TOOLKIT =                SIPS_VERSION 
  ! --------------------------------------------------------------------------

  ! This is the type to store runtime Booleans set and used by the l2cf
  integer, parameter :: RTVSTRINGLENGTH = 1024
  integer, parameter :: RTVARRAYLENGTH  = 128
  
  integer, private :: i ! For loop constructor below

  type :: runTimeValues_T
    ! Two arrays bound as a logical-valued hash
    character(len=RTVSTRINGLENGTH)     :: lkeys = 'true,false'
    logical, dimension(RTVARRAYLENGTH) :: lvalues = &
      & (/ .TRUE., (.FALSE., i=2, RTVARRAYLENGTH) /)
    ! Add two more arrays bound for each kind of hash: integer, string, real, ..
  end type runTimeValues_T
  
  type(runTimeValues_T), save :: runTimeValues
!=============================================================================
contains 
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

END MODULE MLSL2Options
!=============================================================================

!
! $Log$
! Revision 2.35  2006/02/21 19:19:27  pwagner
! New things to create, refer to run time booleans in l2cf
!
! Revision 2.34  2006/02/10 21:13:30  pwagner
! May specify skipRetrivel for particular Phases; dumps may go to special dumpfile
!
! Revision 2.33  2005/07/21 23:40:54  pwagner
! Removed unneeded ILLEGALL1BRADID, MAXNUML1BRADIDS
!
! Revision 2.32  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.31  2005/03/12 00:48:01  pwagner
! Added RESTARTWARNINGS; corrected vsn id
!
! Revision 2.30  2004/12/14 00:04:24  pwagner
! New early stop options added for quicker debugging
!
! Revision 2.29  2004/07/08 22:48:44  pwagner
! Made SIPS_VERSION public
!
! Revision 2.28  2004/04/27 23:49:51  pwagner
! Added SKIPDIRECTWRITES option
!
! Revision 2.27  2004/03/12 00:28:56  pwagner
! At last hdf version at output increased to 5
!
! Revision 2.26  2004/01/23 01:06:39  pwagner
! Added CATENATESPLITS
!
! Revision 2.25  2003/12/05 00:39:35  pwagner
! Added patch option, section timing units
!
! Revision 2.24  2003/11/07 00:46:51  pwagner
! New quicker preflight option: --checkPaths
!
! Revision 2.23  2003/10/09 23:58:34  pwagner
! Updated CURRENT_VERSION_ID to 1.4
!
! Revision 2.22  2003/09/05 23:22:52  pwagner
! Has new SKIPRETRIEVAL option
!
! Revision 2.21  2003/06/09 22:49:32  pwagner
! Reduced everything (PCF, PUNISH.., etc.) to TOOLKIT
!
! Revision 2.20  2003/05/02 20:53:19  pwagner
! Reordered to make SIPS-dependent section clearer; default_hdfversion at read now wildcard
!
! Revision 2.19  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.18  2002/10/03 23:00:03  pwagner
! You can set l1b, l2gp hdfversions on command line
!
! Revision 2.17  2002/08/28 22:25:42  pwagner
! Moved LEVEL1_HDFVERSION, ILLEGALL1BRADID, MAXNUML1BRADIDS here from global_settings
!
! Revision 2.16  2002/03/14 23:38:28  pwagner
! Gets HDFVERSION_4 and 5 from MLSFiles module
!
! Revision 2.15  2002/02/12 00:25:25  pwagner
! New current_version_id parameter
!
! Revision 2.14  2002/02/05 00:44:03  pwagner
! Added garbage collection stuff
!
! Revision 2.13  2002/01/29 23:49:38  pwagner
! Separate DEFAULT_HDFVERSION_(READ)(WRITE)
!
! Revision 2.12  2002/01/23 21:48:16  pwagner
! Added DEFAULT_HDFVERSION
!
! Revision 2.11  2001/09/28 17:57:47  pwagner
! SIPS_VERSION controls other logical options
!
! Revision 2.10  2001/07/16 23:43:15  pwagner
! With settable NORMAL_EXIT_STATUS
!
! Revision 2.9  2001/05/30 22:56:48  pwagner
! Moved PCFL2CFSAMECASE here from OutputAndClose
!
! Revision 2.8  2001/05/15 23:46:08  pwagner
! Removed 2 settings from MLSL2Opts; now in switches
!
! Revision 2.7  2001/05/11 23:48:23  pwagner
! Changed to not echo globals; added note on SIPS
!
! Revision 2.6  2001/05/09 23:34:13  pwagner
! Added ECHO_GLOBAL_STNGS LOG_TO_STDOUT
!
! Revision 2.5  2001/05/06 20:54:40  pwagner
! Default settings should work for most jpl users
!
! Revision 2.4  2001/05/04 22:54:31  pwagner
! Added TOOLKIT, CREATEMETADATA, PCF_FOR_INPUT
!
! Revision 2.3  2001/04/20 20:41:52  pwagner
! Added QUIT_ERROR_THRESHOLD
!
! Revision 2.2  2001/04/17 20:26:28  pwagner
! Added OUTPUT_PRINT_UNIT
!
! Revision 2.1  2001/04/16 23:53:10  pwagner
! First commit
!
