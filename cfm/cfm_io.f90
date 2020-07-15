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
   use machine, only: is_a_directory
   implicit none

   public :: Read_Spectroscopy, ReadDACSFilterShapes, ReadAntennaPatterns
   public :: ReadFilterShapes, ReadPointingGrids, ReadPFAFile, ReadHDF5L2PC
   public :: Read_ptan, Write_To_File1, Write_To_File2
   public :: Write_To_HDF5, Write2HDF5 !added by Zheng Qu
   public :: read_H5inputdata !added by Zheng Qu
   public :: read_txtinputdata, path_join !added by Zheng Qu
   public :: Copy_RadianceJacobian, Filter_RadianceJacobian !added by Zheng Qu
   public :: IDL_where, log_inputdata  !added by Zheng Qu
   public :: runLogFileUnit, errorLogFileUnit, runlog, & ! Added by Zheng Qu
                        errlog, openLogFiles  ! Added by Zheng Qu
   
    !  arguments for  MLSMessageSetup 
    integer , parameter ::  runLogFileUnit =131 ! LogFileUnit is the logical unit for process run log file
    integer , parameter ::  errorLogFileUnit=130 ! LogFileUnit is the logical unit for error log file

   private
   !---------------------------- RCS Ident Info -------------------------------
   character(len=*), parameter :: ModuleName="$RCSfile$"
   !---------------------------------------------------------------------------

   contains
   
       function cstring2fortranchar(strin) result (strout)
       ! This is to convert a C string (ends with null) to Fortran character variable for file name
       ! It replaces all trailing characters since the first occurence of 'null' with blank space
            character(len=*), intent(in) :: strin
            character(len=len(strin)) :: strout
            integer j
            strout = strin
            do j = 1, len(trim(strin))
                if (ichar(strin(j:j)) .eq. 0) then
                    strout(j:len(strin)) = ' '
                endif
            enddo
       
       
       end function cstring2fortranchar
       
       subroutine errlog(mode_name, msgs)
       ! Before calling MLSMessageSetup, the error message doesn't goto specified log file
       ! So we use this subroutine to handle errors before MLSMessageSetup is called.
           use MLSMessageModule, only: MLSMessageExit
           character(len=*), intent(in) :: mode_name
           character(len=*), intent(in), dimension(:) :: msgs
           integer j
           ! 
           write(errorLogFileUnit,'(A)')   'Error module:'
           write(errorLogFileUnit,'(2A)')   '...', mode_name
           write(errorLogFileUnit,'(A)')  'Error:'
           write(errorLogFileUnit,*)  ''
           do j=1, size(msgs)
                write(errorLogFileUnit,'(8x, A)') msgs(j)
           enddo
           call MLSMessageExit(status=1)
       end subroutine errlog
       
       subroutine runlog(mode_name, msg, msg_only)
           
           character(len=*), intent(in) :: mode_name, msg
           logical, intent(in), optional :: msg_only
           ! 
           ! System date time variable
           character(len=8) :: dateStr ! format: 20120930
           character(len=10) :: timeStr ! format: 080000.000
           character(len=5) :: zoneStr  ! format : +0500
           character(len=31) :: now ! output format: '2012-09-30 08:00:00.000 +0500'
           integer :: timeValues(8) ! not used
           
           call DATE_AND_TIME(dateStr, timeStr, zoneStr, timeValues)

           ! format for 'now': '2012-09-30 08:00:00.000 +0500'
           now = dateStr(1:4)//'-'//dateStr(5:6)//'-'//dateStr(7:8) &
            //' '//timeStr(1:2)//':'//timeStr(3:4)//':'//timeStr(5:10) &
            //' ('//zoneStr//')'
            if (present(msg_only)) then
              if (msg_only) then
                write(runLogFileUnit,'(8x, A)')  msg
              else 
                write(runLogFileUnit,'(A,1x,A,": ",A)')  mode_name, now, msg
              endif
            else 
              write(runLogFileUnit,'(A,1x,A,": ",A)')  mode_name, now, msg
            endif
       end subroutine runlog
        
       subroutine  openLogFiles(inputH5fileName, VID, strlen, logProcess, logError)
           use MLSMessageModule, only: MLSMessageExit
           character(len=*) , intent(in) :: inputH5fileName, VID
           integer, intent(in) :: strlen
           character(len=strlen),  intent(out) :: logProcess, logError
           integer :: j1,j0

       
    
        ! Error log file name will be following the same naming pattern as stateVector file name, e.g.
        ! 'TES_0000010310_0768_003_MLS_1681_2009d041_MLSCFM1.57_run.log'
        ! A thorough file name check will be performed after the stateVecor file is read, with runid, profile number etc. available.
           
           j0 = SCAN(inputH5fileName, '/', .TRUE.) ! search for last '/', j0=0 if not found
           j1 = j0 + 1
           
           if (len(inputH5fileName) .lt. 41+j1) then
             print *, 'Preliminary checking: StateVector file name is too short(<42 chars)!'
             print *, 'The file name ' ,inputH5fileName,' did NOT follow the naming convention--'
             print *, 'TES_...MLS_...d*'
             call MLSMessageExit(status=1)
           endif
           
           
           ! construct log file name
           logProcess = adjustl(inputH5fileName(j1:j1+41)//VId//'_run.log')
           logError = adjustl(inputH5fileName(j1:j1+41)//VId//'_error.log')
           print *, 'Process log=',trim(logProcess)
           print *, 'Error log=',trim(logError)
           open(UNIT=errorLogFileUnit,FILE=trim(logError),  ERR=110, STATUS='UNKNOWN')
           open(UNIT=runLogFileUnit,FILE=trim(logProcess), ERR=120, STATUS='UNKNOWN' )
           return 
      
110        print *, 'Unable to open error log file "'//trim(logError)//'" for writing.'
           call MLSMessageExit(status=1)
120        print *, 'Unable to open process run log file "'//trim(logProcess)//'" for writing.'
           call MLSMessageExit(status=1)
        
       end subroutine openLogFiles

       subroutine validate_path(path, modulename, varname)
           character(len=*), intent(in) :: path, modulename, varname
           logical :: existsd, existsf
            
           ! inquire( DIRECTORY=trim(path), exist=existsd )
           inquire( FILE=trim(path), exist=existsf )
           existsd = is_a_directory( path )
           if (.NOT. (existsd .OR. existsf)) then
                call errlog ( moduleName, &
             &(/'Path '//trim(varname)//'='//trim(path)//' is not a valid path!'/)) 
           endif
       end subroutine validate_path
       
       !-------------------
       function path_join(path1, path2, novalidation)
           character(len=*), intent(in) :: path1, path2
           logical, optional, intent(in) :: novalidation
           character(len=len_trim(path1)+len_trim(path2)+1) :: path_join
           integer I, L
           
           character(len(path1))  :: p1
           character(len(path2))  :: p2
           character(len=*), parameter :: moduleName = 'CFM_IO:PATH_JOIN'
           logical :: existsd, existsf
       
           p1 = adjustl(path1)
           p2 = adjustl(path2)
           L = len(trim(p1))
           if (p1(L:L) .eq. '/') then 
               ! path1 contains the trailing '/'
               path_join = trim(p1) // trim(p2) 
           else
               ! path1 does not contain the trailing '/'
               path_join = trim(p1) // '/' // trim(p2) 
           endif
           if (.NOT. present(novalidation) ) then 
               call validate_path(path_join, modulename, 'path_join')
           endif  
       
       end function path_join
    
       function IDL_where(logarr) 
       ! This almost simulates IDL function where: given a logical array logarr, 
       ! it returns an index array with for elements that are true in logarr.
       ! Returning value is a 1-D integer pointer
            use MLSMessageModule, only: MLSMessage, MLSMSG_Error
            integer, dimension(:), pointer :: IDL_where
            logical, dimension(:), intent(in) :: logarr
            integer j, n, nn, nTrue, ierr, i
            character(len=256) CERRMSG
            
           
            nTrue=0
            n = size(logarr)
            do j=1,n
                if(logarr(j)) nTrue=nTrue+1
            enddo
            
            allocate(IDL_where(max(nTrue,1)),STAT=ierr, ERRMSG=CERRMSG)
            if (ierr /=0) then  
                call MLSMessage ( MLSMSG_Error, moduleName, &
            &    'Unable to allocate arrays in IDL_where. ' //CERRMSG)  
                
            endif
            if (nTrue .eq. 0) then
                IDL_where(1) = -1
                return
            endif
            i = 0
            do j=1,n
                if(logarr(j)) then
                  i=i+1
                  IDL_where(i) = j
                endif
            enddo
        
      end function IDL_where   
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
             call Read_HDF5_Spectroscopy(trim(filename))
          else
             call MLSMessage(MLSMSG_Error, moduleName, &
             filetype // " not supported")
          end if
      end subroutine Read_Spectroscopy

!
!

    
! Read some input data from a (text) F90 namelist file, Zheng Qu

    subroutine log_inputdata( TES_MLS_OSP_PATH, MLS_CFM_OUTPUT_PATH, &
        MLS_INPUT_FILE_PATH , fileName, &
        inputH5fileName, spectroscopy, &
        leapsecFile,  antennaPatterns, filterShapes, &
        DACSFilterShapes, pointingGrids, pfaFiles, l2pc , &
            signalFileName, configFileName, &
            vGridStandard37Start ,&
            vGridExtinctionStart ,& 
            hGridStandardVal1 , &
            hGridStandardVal2 , &
            qH2OMinValue , &
            vGridStandard37formula , &
            vGridExtinctionFormula, &
            outputTxtFile, &
            H2O_Sample_Vals, &
            REFGPHINPUT, VGRIDREFGPHVALS, MLS_Year, MLS_Day, &
            MLS_ProfileNumber, Tes_run, Tes_scan, TES_sequnce, &
            l1boa, l1brad, l2GP,l2dgm)
       
       
       use MLSCommon, only: r8, i4
       use MLSMessageModule, only: MLSMessage
       
       character(len=*), intent(in) ::  TES_MLS_OSP_PATH, fileName, &
        inputH5fileName, MLS_CFM_OUTPUT_PATH, MLS_INPUT_FILE_PATH 

     
       character(len=256) ::  spectroscopy, &
        leapsecFile,  antennaPatterns, filterShapes, &
        DACSFilterShapes, pointingGrids, &
        signalFileName, configFileName
       character(len=256), dimension(8) :: pfaFiles
       character(len=256), dimension(1) :: l2pc

       !! Some variables used by mockup with default values
       real(r8) :: vGridStandard37Start 
       real(r8) :: vGridExtinctionStart
       real(r8) :: hGridStandardVal1 
       real(r8) :: hGridStandardVal2 
       real(r8) :: qH2OMinValue 
       character(len=64) :: vGridStandard37formula 
       character(len=64) :: vGridExtinctionFormula
       logical :: outputTxtFile 
       real(r8), dimension(55) :: H2O_Sample_Vals 
        
        !
       real(r8)  REFGPHINPUT, VGRIDREFGPHVALS(1)
       character(len=*) :: L1BOA, L1BRAD, L2GP, L2DGM
      
       integer(i4) :: MLS_Year, MLS_Day, MLS_ProfileNumber, &
        Tes_run,Tes_scan,TES_sequnce
        

       integer :: ierr, j, i
       character (len=255) ::  string_buffer
       character(len=6) :: moduleName='CFM_IO:LOG_INPUTDATA'
       
        ! log input parameters
        
       call runlog ( moduleName, &
        & '***Environment variable for this run***' )

       call runlog ( moduleName, &
        & 'TES_MLS_OSP_PATH='//trim(TES_MLS_OSP_PATH), &
            msg_only=.TRUE.)

       call runlog (  moduleName, &
        & 'MLS_CFM_OUTPUT_PATH='//trim(MLS_CFM_OUTPUT_PATH), &
            msg_only=.TRUE.)
       
       call runlog (  moduleName, &
        & 'MLS_INPUT_FILE_PATH ='//trim(MLS_INPUT_FILE_PATH ), &
            msg_only=.TRUE.)
       
       call runlog (  moduleName, &
        & '***Command line parameter for this run***' )

       call runlog ( moduleName, &
        & 'stateVectorFile='//trim(inputH5fileName) , &
            msg_only=.TRUE.)

       call runlog (moduleName, &
        & '***Parameters read from OSP file ' // &
            trim(path_join(TES_MLS_OSP_PATH,'MLS_CFM/config/'//trim(fileName)))//'***' )
        
       call runlog (  moduleName, &
        & 'spectroscopy='//trim(spectroscopy), &
            msg_only=.TRUE.)
        
       call runlog (  moduleName, &
        & 'leapsecFile='//trim(leapsecFile) , &
            msg_only=.TRUE.)
        
        
       call runlog (  moduleName, &
        & 'antennaPatterns='//trim(antennaPatterns) , &
            msg_only=.TRUE.)
        
       call runlog (  moduleName, &
        & 'filterShapes='//trim(filterShapes) , &
            msg_only=.TRUE.)
        
       call runlog (  moduleName, &
        & 'DACSFilterShapes='//trim(DACSFilterShapes) , &
            msg_only=.TRUE.)
     
       call runlog (  moduleName, &
        & 'pointingGrids='//trim(pointingGrids) , &
            msg_only=.TRUE.)
        
       call runlog (  moduleName, &
        & 'pfaFiles='//trim(pfaFiles(1))//' ,', &
            msg_only=.TRUE.)
       do i = 2, size(pfaFiles)
         call runlog (  moduleName, &
         & trim(pfaFiles(i))//' ,', &
            msg_only=.TRUE.)
       end do

       call runlog (  moduleName, &
        & 'l2pc='//trim(l2pc(1)) , &
            msg_only=.TRUE.)
       do i = 2, size(l2pc)
           call runlog (  moduleName, &
           & trim(l2pc(i))//' ,' , &
            msg_only=.TRUE.)
       end do
        
       call runlog (  moduleName, &
        & 'signalFileName='//trim(signalFileName) , &
            msg_only=.TRUE.)
        
       call runlog (  moduleName, &
        & 'configFileName='//trim(configFileName) , &
            msg_only=.TRUE.)
        
       WRITE(string_buffer,*) vGridStandard37Start
       call runlog (  moduleName, &
        & 'vGridStandard37Start='//trim(string_buffer) , &
            msg_only=.TRUE.)
        
       WRITE(string_buffer,*) vGridExtinctionStart
       call runlog (  moduleName, &
        & 'vGridExtinctionStart='//trim(string_buffer) , &
            msg_only=.TRUE.)
           
       WRITE(string_buffer,*) hGridStandardVal1
       call runlog (  moduleName, &
        & 'hGridStandardVal1='//trim(string_buffer) , &
            msg_only=.TRUE.)
        
        
       WRITE(string_buffer,*) hGridStandardVal2
       call runlog (  moduleName, &
        & 'hGridStandardVal2='//trim(string_buffer), &
            msg_only=.TRUE.)
        
       call runlog (  moduleName, &
        & 'vGridStandard37formula='//trim(vGridStandard37formula) , &
            msg_only=.TRUE.)
        
       call runlog (  moduleName, &
        & 'vGridExtinctionFormula='//trim(vGridExtinctionFormula) , &
            msg_only=.TRUE.)
        
       WRITE(string_buffer,*) outputTxtFile
       call runlog (  moduleName, &
        & 'outputTxtFile='//trim(string_buffer) , &
            msg_only=.TRUE.)
        
        
       call runlog (  moduleName, &
        & 'H2O_Sample_Vals=', &
            msg_only=.TRUE.)
        
        
       do i = 1, size(H2O_Sample_Vals),4
           WRITE(string_buffer,'(4E15.5)') (H2O_Sample_Vals(j), &
               j=i,min(i+3,size(H2O_Sample_Vals)))
           call runlog (  moduleName, &
        &  trim(string_buffer), &
            msg_only=.TRUE.)
       end do
        
                
       call runlog ( moduleName, &
        & '***Parameters read from stateVector file ' // &
            inputH5fileName//'***' )
        
       call runlog (  moduleName, &
        & 'l1boa='//trim(l1boa) , &
            msg_only=.TRUE.)
        
       call runlog (  moduleName, &
        & 'l1brad='//trim(l1brad) , &
            msg_only=.TRUE.)
        
       call runlog (  moduleName, &
        & 'l2GP='//trim(l2GP) , &
            msg_only=.TRUE.)
        
       call runlog (  moduleName, &
        & 'l2dgm='//trim(l2dgm), &
            msg_only=.TRUE.)
        

       WRITE(string_buffer,*) REFGPHINPUT
       call runlog (  moduleName, &
        & 'REFGPHINPUT='//trim(string_buffer) , &
            msg_only=.TRUE.)
        
        
       WRITE(string_buffer,*) VGRIDREFGPHVALS
       call runlog (  moduleName, &
        & 'VGRIDREFGPHVALS='//trim(string_buffer) , &
            msg_only=.TRUE.)
        
        
       WRITE(string_buffer,*) MLS_ProfileNumber
       call runlog (  moduleName, &
        & 'MLS_ProfileNumber='//trim(string_buffer) , &
            msg_only=.TRUE.)
        
        
       WRITE(string_buffer,*) MLS_Year
       call runlog (  moduleName, &
        & 'MLS_Year='//trim(string_buffer) , &
            msg_only=.TRUE.)
        
       WRITE(string_buffer,*) Tes_run
       call runlog (  moduleName, &
        & 'Tes_run='//trim(string_buffer), &
            msg_only=.TRUE.)
        
        
       WRITE(string_buffer,*) Tes_scan
       call runlog (  moduleName, &
        & 'Tes_scan='//trim(string_buffer) , &
            msg_only=.TRUE.)
        
        
       WRITE(string_buffer,*) TES_sequnce
       call runlog (  moduleName, &
        & 'TES_sequnce='//trim(string_buffer) , &
            msg_only=.TRUE.)
            
    end subroutine log_inputdata
  

! Read some input data from a (text) F90 namelist file, Zheng Qu
    subroutine read_txtinputdata(MLS_CFM_OSP_PATH, fileName, spectroscopy, &
        leapsecFile,  antennaPatterns, filterShapes, &
        DACSFilterShapes, pointingGrids, pfaFiles, l2pc , &
            signalFileName, configFileName, &
            vGridStandard37Start ,&
            vGridExtinctionStart ,& 
            hGridStandardVal1 , &
            hGridStandardVal2 , &
            qH2OMinValue , &
            vGridStandard37formula , &
            vGridExtinctionFormula, &
            outputTxtFile, &
            H2O_Sample_Vals)
       
       
       use MLSCommon, only: r8

       

       character(len=*), intent(in) ::  MLS_CFM_OSP_PATH, fileName

       integer :: ierr, j
       character (len=255) ::  CERRMSG
     
       character(len=256) ::  spectroscopy, &
        leapsecFile,  antennaPatterns, filterShapes, &
        DACSFilterShapes, pointingGrids, &
        signalFileName, configFileName
       character(len=256), dimension(8) :: pfaFiles
       character(len=256), dimension(1) :: l2pc

       !! Some variables used by mockup with default values
        real(r8) :: vGridStandard37Start 
        real(r8) :: vGridExtinctionStart
        real(r8) :: hGridStandardVal1 
        real(r8) :: hGridStandardVal2 
        real(r8) :: qH2OMinValue 
        character(len=64) :: vGridStandard37formula 
        character(len=64) :: vGridExtinctionFormula
        logical :: outputTxtFile 
        real(r8), dimension(55) :: H2O_Sample_Vals 
        character(len=len(MLS_CFM_OSP_PATH)+7) MLS_CFM_OSP_PATH_config
        character(len=len(MLS_CFM_OSP_PATH)+6) MLS_CFM_OSP_PATH_L2Cal
       
       
    ! Use the f90 standard namelist to initialize parameters
    !
    ! Note that NAMELIST syntax is similar to COMMON BLOCKs
    
        namelist /STATIC_VALS/ vGridStandard37Start ,&
            vGridExtinctionStart ,& 
            hGridStandardVal1 , &
            hGridStandardVal2 , &
            qH2OMinValue , &
            vGridStandard37formula , &
            vGridExtinctionFormula, &
            outputTxtFile
            
        namelist /file_names1/ spectroscopy, &
        leapsecFile,  antennaPatterns, filterShapes, &
        DACSFilterShapes, pointingGrids, &
        signalFileName, configFileName
        
        namelist /file_names2/ pfaFiles
        namelist /file_names3/ l2pc
        namelist /H2O_VALS/ H2O_Sample_Vals
        
        

    ! the OPEN statement defines many of the
    ! NAMELIST characteristics
        ! need the delim, else some implementations will not surround
        ! character strings with delimiters
        ! recl limits the I/O to 160 character lines
        MLS_CFM_OSP_PATH_config = path_join(MLS_CFM_OSP_PATH,'config')
        MLS_CFM_OSP_PATH_L2Cal = path_join(MLS_CFM_OSP_PATH,'L2Cal')
        open(8,file=trim(path_join(MLS_CFM_OSP_PATH_config, &
            trim(adjustl(fileName)))), status='OLD', recl=160, delim='APOSTROPHE')
        read(8,nml=STATIC_VALS)

        read(8,nml=file_names1)
        spectroscopy = path_join(MLS_CFM_OSP_PATH_L2Cal, spectroscopy)
        leapsecFile = path_join('./in', leapsecFile)
        antennaPatterns = path_join(MLS_CFM_OSP_PATH_L2Cal, antennaPatterns)
        filterShapes = path_join(MLS_CFM_OSP_PATH_L2Cal, filterShapes)
        DACSFilterShapes = path_join(MLS_CFM_OSP_PATH_L2Cal, DACSFilterShapes)
        pointingGrids = path_join(MLS_CFM_OSP_PATH_L2Cal, pointingGrids)

        signalFileName = path_join(MLS_CFM_OSP_PATH_config, signalFileName)
        configFileName = path_join(MLS_CFM_OSP_PATH_config, configFileName)
        
        read(8,nml=file_names2)
        do j=1, size(pfaFiles)
            pfaFiles(j) = path_join(MLS_CFM_OSP_PATH_L2Cal, pfaFiles(j))
        enddo
        
        read(8,nml=file_names3)
        do j=1, size(l2pc)
            l2pc(j) = path_join(MLS_CFM_OSP_PATH_L2Cal, l2pc(j))
        enddo
        
        read(8,nml=H2O_VALS)
        

        close(8)
        
    end subroutine read_txtinputdata
    
    
    
    subroutine copy_float32Tofloat64(ptr_from, ptr_to, varname, all_positive)
        ! Copy from a 32bit float array to 64bit array
        ! Both input and output are pointers.
        ! Input pointer will be nullified upon return
        use MLSMessageModule, only: MLSMessage, MLSMSG_Error
        use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
        use MLSCommon, only: r8, r4, i4
        real(r4), dimension(:), pointer :: ptr_from ! ptr_from will be nullified upon return
        real(r8), dimension(:), pointer :: ptr_to
        character(len=*), intent(in) :: varname
        logical, optional :: all_positive
        
        ! Local variables
        character(len=128) :: moduleName='CFM_IO:copy_float32Tofloat64', &
            CERRMSG
        integer :: iostat, ierr, nn, j
        integer, dimension(:), pointer ::  ind_positive=>NULL()

      !----------------
      

        nn = size(ptr_from)
        
        if (present(all_positive)) then 
            ind_positive => IDL_where(ptr_from .GT. 0.)
      
            if (ind_positive(1) .lt. 0) then 
                call MLSMessage ( MLSMSG_Error, moduleName, &
                & 'All elements in '//trim(varname)//' are 0 or negative. MLSCFM failed.' )
            endif
            
            nn = size(ind_positive)

            allocate(ptr_to(nn), STAT=ierr, ERRMSG=CERRMSG)
            if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
                    & 'Unable to allocate array for '//trim(varname)//&
                    ' ; '//CERRMSG)
                
            do j=1,nn
                ptr_to(j) = ptr_from(ind_positive(j))
            enddo
            call Deallocate_test ( ind_positive, &
                & 'Unable to deallocate array for ind_positive when copying '//trim(varname), ModuleName )

        else

            allocate(ptr_to(nn), STAT=ierr, ERRMSG=CERRMSG)
            if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
                    & 'Unable to allocate array for '//trim(varname)//&
                    ' ; '//CERRMSG)
                
            ptr_to = ptr_from
        endif
        call Deallocate_test ( ptr_from, &
                & 'Unable to deallocate array for '//trim(varname), ModuleName )

        return 
    end subroutine copy_float32Tofloat64
      
    subroutine read_H5inputdata (MLS_INPUT_FILE_PATH , fileName, &
        TemperatureInput, H2OInput , &
        O3Input, SO2Input, HNO3Input, COInput, extinctionV2R3Input, &
        PressureCOInput, PressureStandardInput, &
        REFGPHINPUT, VGRIDREFGPHVALS, MLS_Year, MLS_Day,&
        MLS_ProfileNumber, Tes_run,Tes_scan,TES_sequnce, &
        l1boa, l1brad, l2GP,l2dgm, H2O_Sample_Vals)
      ! Read stateVector files (hdf5) for MLS atmospheric profiles, TES apriori and MLS file names, as well as 
      ! TES run-seq-scan + MLS year and day of year etc.
      ! 
      ! Please note the stateVector files are written with values in float32
      ! while MLSCFM calling inputs are float64 (double precision)
      ! We need to convert the arrays before returning
      ! 
      ! MLS file base names are read from stateVector file, we need to add MLS_INPUT_FILE_PATH  to construct full paths
      ! Also need to identify which is which (L1BOA, L1BRAD, L2GP, and L2DGM) for the four MLS files.
      !
      use MLSCommon, only: r8, r4, i4
      use MLSHDF5, only: loadPtrFromHDF5DS, LoadFromHDF5DS
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error
      use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
      use HDF5, only:  H5FOpen_F, H5FClose_F, HSize_T, &
            H5F_ACC_RDONLY_F, h5open_f
      
      
      character(len=*), intent(in) :: fileName, MLS_INPUT_FILE_PATH 
      real(r8), dimension(55), intent(in) :: H2O_sample_vals 

      real(r8), dimension(:), pointer ::  TemperatureInput, H2OInput , &
        O3Input, SO2Input, HNO3Input, COInput, extinctionV2R3Input, &
        PressureCOInput, PressureStandardInput
      real(r8)  REFGPHINPUT, VGRIDREFGPHVALS(1)
      character(len=*) :: L1BOA, L1BRAD, L2GP, L2DGM
      
      real(r4), dimension(:), pointer ::  TemperatureInput32, H2OInput32 , &
        O3Input32, SO2Input32, HNO3Input32, COInput32, extinctionV2R3Input32, &
        PressureCOInput32, PressureStandardInput32, &
        REFGPHINPUT32 , VGRIDREFGPHVALS32
      integer(i4), dimension(1) :: MLS_Year32, MLS_Day32, &
        MLS_ProfileNumber32, &
        Tes_run32, Tes_scan32, TES_sequence32
      
      integer(i4) :: MLS_Year, MLS_Day, MLS_ProfileNumber, &
        Tes_run,Tes_scan,TES_sequnce
      
        
      ! Local variables
      character(len=128) :: moduleName="CFM_IO:read_H5inputdata", CERRMSG

      integer :: iostat, fileID, error, groupID, nH2OInput, nn, ierr, j
      integer, parameter :: strlen=512, nfiles=5
      character(len=strlen) :: MLS_Files(nfiles)
      character(len=6) :: MLS_Files_missing(4) = &
                (/'L1BOA ', 'L1BRAD', 'L2GP  ', 'L2DGM '/)
      character(len=10) :: YYYYDDD 
      

      ! Open HDF file
      call h5fopen_f ( trim(fileName), H5F_ACC_RDONLY_F, fileID, iostat )
      if (iostat /= 0) call MLSMessage (MLSMSG_Error, ModuleName, &
         & 'Unable to open input HDF5 file ' // trim(fileName) // '.' )


      call loadPtrFromHDF5DS ( fileID, '/Data/MLS_TATM_Values', TemperatureInput32 )
      !call loadPtrFromHDF5DS ( fileID, 'H2OInput', H2OInput )
      call loadPtrFromHDF5DS ( fileID, '/Data/MLS_O3_Values', O3Input32 )
      call loadPtrFromHDF5DS ( fileID, '/Data/MLS_SO2_Values', SO2Input32 )
      call loadPtrFromHDF5DS ( fileID, '/Data/MLS_HNO3_Values', HNO3Input32 )


      call loadPtrFromHDF5DS ( fileID, '/Data/TES_CO_Apriori', COInput32 )
      call loadPtrFromHDF5DS ( fileID, '/Data/MLS_Extinction_Values', extinctionV2R3Input32 )
      call loadPtrFromHDF5DS ( fileID, '/Data/TES_Pressure_Values', PressureCOInput32 )
      call loadPtrFromHDF5DS ( fileID, '/Data/MLS_TATM_Pressure_Values', PressureStandardInput32 )

      call loadPtrFromHDF5DS ( fileID, '/Data/MLS_refGPH_Values', REFGPHINPUT32 )
      call loadPtrFromHDF5DS ( fileID, '/Data/MLS_refGPH_Pressure_Values', VGRIDREFGPHVALS32 )

      call LoadFromHDF5DS ( fileID, '/Matched/MLS_Year', MLS_Year32 )
      call LoadFromHDF5DS ( fileID, '/Matched/MLS_Day', MLS_Day32 )

      call LoadFromHDF5DS ( fileID, '/Matched/MLS_ProfileNumber', MLS_ProfileNumber32 )
      call LoadFromHDF5DS ( fileID, '/Matched/TES_Run', TES_Run32 )
      call LoadFromHDF5DS ( fileID, '/Matched/TES_Scan', TES_Scan32 )
      call LoadFromHDF5DS ( fileID, '/Matched/TES_Sequence', TES_Sequence32 )
      
      call LoadFromHDF5DS (fileID, '/Matched/MLS_Files', MLS_Files)
      

      ! string for MLS year and day number in file name. e.g., '_2009d041_'
      write(YYYYDDD, '("_",I4.4,"d",I3.3,"_")') MLS_Year32(1), MLS_Day32(1)
      
      !
      ! MLS_Files should contain at least 4 file names for L1BOA, L1BRAD, L2GP, and L2DGM
      ! We need to identify which is which and assign the full path to corresponding vars.
      !

      do j=1, nfiles
        if (INDEX(MLS_Files(j), YYYYDDD) .LE. 0 ) then
            continue ! year day number not match.
        endif
        
        if (INDEX(MLS_Files(j), '_L1BOA_') .GT. 0) then
            L1BOA = path_join(MLS_INPUT_FILE_PATH , &
                trim(adjustl(cstring2fortranchar(MLS_Files(j)))))
            MLS_Files_missing(1) = ' '
        else if (INDEX(MLS_Files(j), '_L1BRADG_') .GT. 0) then
            L1BRAD = path_join(MLS_INPUT_FILE_PATH , &
                trim(adjustl(cstring2fortranchar(MLS_Files(j)))))
            MLS_Files_missing(2) = ' '
        else if (INDEX(MLS_Files(j), '_L2GP-CO_') .GT. 0) then
            L2GP = path_join(MLS_INPUT_FILE_PATH , &
                trim(adjustl(cstring2fortranchar(MLS_Files(j)))))
            MLS_Files_missing(3) = ' '
        else if (INDEX(MLS_Files(j), '_L2AUX-DGM_') .GT. 0) then
            L2DGM = path_join(MLS_INPUT_FILE_PATH , &
                trim(adjustl(cstring2fortranchar(MLS_Files(j)))))
            MLS_Files_missing(4) = ' '
        endif
      enddo

      !In case one MLS file is not found in  MLS_Files data field.
      do j = 1, 4
        if (len_trim(MLS_Files_missing(j)) .NE. 0 ) then

           call errlog(ModuleName,  &
           & (/'Unable to locate ' // trim(MLS_Files_missing(j))// &
           ' file from "/Matched/MLS_Files" dataset in stateVector file' , &
           trim(fileName) // '!' /))
        endif
      enddo 
      
      ! check if path is valid 
      call validate_path(L1BOA, moduleName, 'L1BOA')
      call validate_path(L1BRAD, moduleName, 'L1BRAD')
      call validate_path(L2GP, moduleName, 'L2GP')
      call validate_path(L2DGM, moduleName, 'L2DGM')
      
      ! convert to float64
      
      call copy_float32Tofloat64 ( TemperatureInput32, TemperatureInput, 'TemperatureInput' )
      call copy_float32Tofloat64 ( O3Input32,O3Input, 'O3Input')
      call copy_float32Tofloat64 ( SO2Input32, SO2Input, 'SO2Input' )
      call copy_float32Tofloat64 ( HNO3Input32, HNO3Input, 'HNO3Input' )

      call copy_float32Tofloat64 ( COInput32, COInput, 'COInput' , all_positive=.TRUE.)

      call copy_float32Tofloat64 ( extinctionV2R3Input32, extinctionV2R3Input, 'extinctionV2R3Input' )

      call copy_float32Tofloat64 ( PressureCOInput32, PressureCOInput, 'PressureCOInput',  all_positive=.TRUE. )

      
      call copy_float32Tofloat64 ( PressureStandardInput32, PressureStandardInput, 'PressureStandardInput' )
      
      REFGPHINPUT = REFGPHINPUT32(1)
      MLS_Year = MLS_Year32(1)
      MLS_Day = MLS_Day32(1)
      MLS_ProfileNumber = MLS_ProfileNumber32(1)
      VGRIDREFGPHVALS = VGRIDREFGPHVALS32(1)
      Tes_run = Tes_run32(1)
      Tes_scan = Tes_scan32(1)
      TES_sequnce = TES_sequence32(1)
      
      
      call Deallocate_test ( REFGPHINPUT32, &
                & 'Unable to deallocate array for REFGPHINPUT32', ModuleName )

      
      
      call Deallocate_test ( VGRIDREFGPHVALS32, &
                & 'Unable to deallocate array for VGRIDREFGPHVALS32', ModuleName )


      
      nH2OInput = size(TemperatureInput) ! H2O use the same dimension as TATM
      allocate(H2OInput(nH2OInput),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for H2OInput.' )
      
      ! fill H2OInput with sample values (read from osp file in H2O_sample_vals) 
      ! -- it is needed for calling CFm but will actually have no effects 
      ! since PTAN values read form MLS DGG file will ovverride that computed using H2O.
      nn =  min(size(H2O_sample_vals), nH2OInput)
      H2OInput(1:nn) = H2O_sample_vals(1:nn)
      if (nH2OInput .GT. nn) then
          H2OInput(nn+1:nH2OInput) = H2OInput(nn) ! use a small value to fill up
      endif
      

      ! close file
      call H5FClose_F ( fileID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to close hdf5 input file ' // trim(fileName) // '.' )   

   end subroutine  read_H5inputdata

   ! Read ptan values from L2AUX-DGM file, Pranjit Saha
   subroutine Read_ptan ( filename, fileType, MAFnumber, ptanValuesRead, &
     & DSName )
      use MLSHDF5, only: LoadFromHDF5DS, LoadPtrFromHDF5DS, GetHDF5DSDims, &
                         IsHDF5DSPresent
      use MLSCommon, only: r8
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error, &
                                  MLSMSG_Allocate, MLSMSG_DeAllocate
      use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
      use Parse_Signal_m, only: Parse_Signal
      use HDF5, only: H5F_ACC_RDONLY_F, H5FOpen_F, H5FClose_F, HSize_T

      character(len=*), intent(in)           :: filename
      character(len=*), intent(in)           :: filetype
      integer, intent(in)                    :: MAFnumber
      real(r8), dimension(125)               :: ptanValuesRead
      character(len=*), optional, intent(in) :: DSName

      ! Local variables
      integer :: iostat, fileID, Line1, nLines, lineN, j, i
      integer(hsize_t) :: Shp(1), Shp2(2) ! To get the shapes of datasets HD
      real(r8), pointer :: ptanValues(:,:,:)
      character(len=128) :: PTanName

      ! Open HDF file
      nullify( ptanValues )
      PTanName = 'GHz.ptan-FinalPtan'
      if ( present(DSName) ) PTanName = DSName
      call h5fopen_f ( trim(fileName), H5F_ACC_RDONLY_F, fileID, iostat )
      if (iostat /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
         & 'Unable to open HDF5 L2AUX-DGM file ' // trim(fileName) // '.' )

      ! Check if field (GHz.ptan-FinalPtan) is present in the file
      if ( IsHDF5DSPresent ( fileID, trim(PTanName) ) ) then
        ! Read data from HDF file
        call  loadPtrFromHDF5DS ( fileID, trim(PTanName), ptanValues )
        do j = 1, 125
           ptanValuesRead(j) = ptanValues(1, j, MAFnumber)
        end do
      end if

      ! clode file and deallocate memory
      call H5FClose_F ( fileID, iostat )

   end subroutine Read_ptan

   ! Write 'band9', 'precision9' in separate text files, Pranjit Saha
   subroutine Write_To_File2 ( band9, precision9 )
      use VectorsModule, only: VectorValue_T, DumpMask
      use VectorsModule, only: VECTOR_T, VectorValue_T
      use QuantityTemplates, only: QUANTITYTEMPLATE_T
      use MLSMessageModule, only:  MLSMessage, MLSMSG_Error

      type(VectorValue_T), intent(in) :: band9
      type(VectorValue_T), intent(in) :: precision9

      integer :: columns
      integer :: rows
      integer :: indx, i
      CHARACTER(LEN=20), PARAMETER :: FMT1 = "(F16.5)"
      CHARACTER(LEN=20), PARAMETER :: FMT2 = "(E18.8)"

      character(len=*), parameter :: ModuleName = 'CFM_IO:Write_To_File2'
      
      ! Open file to write band9
      open (unit = 1, file = "band9.txt")
      write (1, *) "Number of columns: ", band9%template%noChans
      write (1, *) "Number of rows: ",    band9%template%noSurfs
      do rows = 1, band9%template%noChans * band9%template%noSurfs
         write (1, FMT1, advance="no") band9%values(rows, 1)
         ! Advance a line after writing noSurfs values
         i = mod(rows, band9%template%noChans)
         if ( i == 0 ) then
             write (1, *) ""
         endif
      end do
      close(1)

      ! Open file to write precision9
      open (unit = 1, file = "precision9.txt")
      write (1, *) "Number of columns: ", precision9%template%noChans
      write (1, *) "Number of rows: ", precision9%template%noSurfs
      do rows = 1, precision9%template%noChans * precision9%template%noSurfs
         write (1, FMT1, advance="no") precision9%values(rows, 1)
         ! Advance a line after writing noChans values
         i = mod(rows, precision9%template%noChans)
         if ( i == 0 ) then
            write (1, *) ""
         endif
      end do
      close(1)

      ! Code from Paul Wagner, See mail on Feb 24, 2012
      ! Add DumpMask to "use VectorsModule, only: .."  statement
      if ( .not. associated(precision9%mask) ) then
         call MLSMessage (MLSMSG_Error, ModuleName, &
         & 'Sorry--mask not associated in precision9.' )
      endif
      !Open file to write mask
      open (unit = 1, file = "mask9.txt")
      write (1, *) "Number of columns: ", precision9%template%noChans
      write (1, *) "Number of rows: ", precision9%template%noSurfs
      do rows = 1, precision9%template%noChans * precision9%template%noSurfs
          !! write (1, format='(z3.2)', advance="no") precision9%mask(rows, 1)

            write (1, fmt='(z3.2)', advance="no") ichar(precision9%mask(rows, 1))

          ! Advance a line after writing noChans values
          i = mod(rows, precision9%template%noChans)
          if ( i == 0 ) then
             write (1, *) ""
          endif
      end do
      close(1)

      !call DumpMask( precision9, details=0 )
      !call DumpMask( precision9, details=2 )

   end subroutine  !! end of subroutine Write_To_File2 ( band9, precision9 )

   ! Write 'state;, 'radiance', 'jacobian', 'ptanGHz' to separate text file, Pranjit Saha
   ! You will be happier if you rewrite this subroutine so it
   ! doesn't require so much tweaking each time you add or remove
   ! a species from state
   subroutine Write_To_File1 (state, radiance, jacobian, newer_ptanGHz)
      ! Reference to user defined types
      use VectorsModule, only: VECTOR_T, VectorValue_T
      use MATRIXMODULE_1 , only: MATRIX_T

      ! Function argument type
      type(Vector_T),           intent(in) :: state
      type(Vector_T),           intent(in) :: radiance
      type(Matrix_T),           intent(in) :: jacobian
      type(VectorValue_T), intent(in) :: newer_ptanGHz

      ! Local variables
      integer :: columns
      integer :: rows
      integer :: indx, i
      CHARACTER(LEN=20), PARAMETER :: FMT1 = "(F16.5)"
      CHARACTER(LEN=20), PARAMETER :: FMT2 = "(E18.8)"
      integer :: qty

      ! Open file to write, for 'state'
      qty = 0
      open (unit = 1, file = "state.txt")
      !write (1, *) ""
      qty = qty + 1
      write (1, *) "Species 1: Temperature"
      write (1, *) "Total values: ", state%quantities(qty)%template%noSurfs

      do rows = 1, state%quantities(qty)%template%noSurfs
        write (1, FMT1, advance="no") state%quantities(qty)%values(rows, 1)
      end do

      !write (1, *) ""
      ! You omitted CO, remember
      if ( .false. ) then
      qty = qty + 1
      write (1, *) "Species 2: CO"
      write (1, *) "Total values: ", state%quantities(qty)%template%noSurfs
      do rows = 1, state%quantities(qty)%template%noSurfs
        write (1, FMT2, advance="no") state%quantities(qty)%values(rows, 1)
      end do
      endif

      !write (1, *) ""
      qty = qty + 1
      write (1, *) "Species 3: O2"
      write (1, *) "Total values: ", state%quantities(qty)%template%noSurfs
      do rows = 1, state%quantities(qty)%template%noSurfs
         write (1, FMT1, advance="no") state%quantities(qty)%values(rows, 1)
      end do

      !write (1, *) ""
      qty = qty + 1
      write (1, *) "Species 4: SO2"
      write (1, *) "Total values: ", state%quantities(qty)%template%noSurfs
      do rows = 1, state%quantities(qty)%template%noSurfs
         write (1, FMT2, advance="no") state%quantities(qty)%values(rows, 1)
      end do

      !write (1, *) ""
      qty = qty + 1
      write (1, *) "Species 5: HNO3"
      write (1, *) "Total values: ", state%quantities(qty)%template%noSurfs
      do rows = 1, state%quantities(qty)%template%noSurfs
         write (1, FMT2, advance="no") state%quantities(qty)%values(rows, 1)
      end do

      !write (1, *) ""
      qty = qty + 1
      write (1, *) "Species 6: O3"
      write (1, *) "Total values: ", state%quantities(qty)%template%noSurfs
      do rows = 1, state%quantities(qty)%template%noSurfs
         write (1, FMT2, advance="no") state%quantities(qty)%values(rows, 1)
      end do

      !write (1, *) ""
      qty = qty + 1
      write (1, *) "Species 7: ExtinctionV2R3"
      write (1, *) "Total values: ", state%quantities(qty)%template%noSurfs
      do rows = 1, state%quantities(qty)%template%noSurfs
         write (1, FMT2, advance="no") state%quantities(qty)%values(rows, 1)
      end do
      ! close file opened to write
      close(1)

      ! To write 'radiance' values in a 2 dimensional format
      ! Dimension: rows = noSurfs, columns = noChans
      ! Open file to write
      open (unit = 1, file = "simulatedRadiance.txt")
      write (1, *) "Number of rows: ", radiance%quantities%template%noChans
      write (1, *) "Number of columns: ", radiance%quantities%template%noSurfs

      do rows = 1, radiance%quantities(1)%template%noChans * radiance%quantities(1)%template%noSurfs
         write (1, FMT1, advance="no") radiance%quantities(1)%values(rows, 1)
         ! Advance a line after writing noChans values
         i = mod(rows, radiance%quantities(1)%template%noChans)
         if ( i == 0 ) then
            write (1, *) ""
         endif
      end do
      ! close file opened to write
      close(1)

      ! Save 'jacobian' values in a file
      open (unit = 1, file = "jacobian.txt")

      rows = 1
      ! Don't you remember? There is no CO
      if ( .false. ) then
      write (1, *) "Species 2: CO"
      write (1, *) "Rows: ", jacobian%block(1, 2)%nRows
      write (1, *) "Columns: ", jacobian%block(1, 2)%nCols

      columns = 2   ! CO is saved in index 2, temp is in index 1 - if retrieved
      do i = 1, jacobian%block(1, 2)%nRows
         do indx = 1, jacobian%block(1, 2)%nCols
           write (1, FMT1, advance="no") jacobian%block(rows, columns) %values(i, indx)
         end do
         write(1, *) ""
      end do
      endif

      columns = 5 ! 6 ! for O3
      write (1, *) "Species 6: O3"
      write (1, *) "Rows: ", jacobian%block(1, columns)%nRows
      write (1, *) "Columns: ", jacobian%block(1, columns)%nCols
      do i = 1, jacobian%block(1, columns)%nRows
         do indx = 1, jacobian%block(1, columns)%nCols
           write (1, FMT1, advance="no") jacobian%block(rows, columns)%values(i, indx)
         end do
         write(1, *) ""
      end do

      columns = 6 ! 7 ! for EXTINCTIONV2
      write (1, *) "Species 7: ExtinctionV2"
      write (1, *) "Rows: ", jacobian%block(1, columns)%nRows
      write (1, *) "Columns: ", jacobian%block(1, columns)%nCols
      do i = 1, jacobian%block(1, columns)%nRows
         do indx = 1, jacobian%block(1, columns)%nCols
           write (1, FMT1, advance="no") jacobian%block(rows, columns)%values(i, indx)
         end do
         write(1, *) ""
      end do
      ! close file opened to write
      close(1)

      ! Write ptanGHz in a text file
      ! Open file to write ptanGHz
      open (unit = 1, file = "ptanGHz.txt")
      write (1, *) "Number of columns: ",  newer_ptanGHz%template%noSurfs
      write (1, *) "Number of rows: ",     newer_ptanGHz%template%noChans
      do rows = 1, newer_ptanGHz%template%noChans * newer_ptanGHz%template%noSurfs
         write (1, FMT1, advance="no") newer_ptanGHz%values(rows, 1)
      end do
      close(1)

   end subroutine     !end of subroutine Write_To_File1(...)
!

   ! Write 'state;, 'radiance', 'jacobian', 'ptanGHz' to a HDF5 file, Zheng Qu
   subroutine Copy_RadianceJacobian ( radiance, jacobian, newer_ptanGHz, &
        Radiance_Calculated, Jacobian_Concatenated, ptanGHz_save)
      ! Reference to user defined types
      use VectorsModule, only: VECTOR_T, VectorValue_T
      use MATRIXMODULE_1 , only: MATRIX_T
      use MLSKINDS, only: RV
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error

      ! Function argument type
      type(Vector_T),           intent(in) :: radiance
      type(Matrix_T),           intent(in) :: jacobian

      type(VectorValue_T), intent(in) :: newer_ptanGHz
      real(rv), dimension(:,:), pointer :: &
            Radiance_Calculated, &
            Jacobian_Concatenated , &
            ptanGHz_save
            
      ! Local variables
      character(len=128) :: moduleName="CFM_IO:Copy_RadianceJacobian"
      integer :: iostat, fileID, error, groupID, ierr
      
      ! Local variables
      integer :: columns, ncolumns
      integer :: rows, nrows
      integer :: indx, i
      integer :: nRows_CO, ncolumns_CO, nRows_O3, ncolumns_O3, &
        nRows_Ext, ncolumns_Ext, ncolumns_all
      
      character(len=128) :: CERRMSG
      
      
      !   get ptanGHZ
      
      !write (1, *) "Number of columns: ",  newer_ptanGHz%template%noSurfs
      !write (1, *) "Number of rows: ",     newer_ptanGHz%template%noChans
      ! Note, in old version, when writing ptanGHz.txt file, # of columns and rows are switched
      nrows = newer_ptanGHz%template%noSurfs
      ncolumns  = newer_ptanGHz%template%noChans
      allocate(ptanGHz_save(nrows, ncolumns),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for ptanGHz_save.'//CERRMSG )
      ptanGHz_save = transpose(reshape(newer_ptanGHz%values(1:nrows*ncolumns, 1), &
        (/ncolumns, nrows/)))
      ! ptanGHz dimension = (1,125)
      
      !open (unit = 1, file = "simulatedRadiance.txt")
      !write (1, *) "Number of rows: ", radiance%quantities%template%noChans
      !write (1, *) "Number of columns: ", radiance%quantities%template%noSurfs
      ! Note, in old version, when writing simulatedRadiance.txt file, # of columns and rows are switched
      
      ncolumns  = radiance%quantities(1)%template%noChans
      nrows = radiance%quantities(1)%template%noSurfs

      allocate(Radiance_Calculated(nrows, ncolumns),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for Radiance_Calculated.'//CERRMSG )
      Radiance_Calculated = &
        transpose(RESHAPE( radiance%quantities(1)%values(1:nrows*ncolumns, 1), &
        (/ ncolumns , nrows/) ))
        

      ! Concatenate 'Jacobian' matrices 
      ! CO_JACOBIAN     DOUBLE    = Array[3125, 65]
      nRows_CO = jacobian%block(1, 2)%nRows
      ncolumns_CO = jacobian%block(1, 2)%nCols
      
      nRows_O3 = jacobian%block(1, 6)%nRows
      ncolumns_O3 = jacobian%block(1, 6)%nCols
      
      nRows_Ext = jacobian%block(1, 7)%nRows
      ncolumns_Ext = jacobian%block(1, 7)%nCols
      
      ncolumns_all = ncolumns_CO+ncolumns_O3+ncolumns_Ext
      
    ! jacob_all = fltarr(nnCO[0],nnCO[1]+nnO3[1]+nnExt[1])
      allocate(Jacobian_Concatenated(nRows_CO, &
        ncolumns_all),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
       & 'Unable to allocate arrays for Jacobian_Concatenated.'//CERRMSG)
        
    !jacob_all[*,0:nnCO[1]-1] = CO_JACOBIAN
       Jacobian_Concatenated(1:nRows_CO,1:ncolumns_CO) = &  
            jacobian%block(1, 2)%values(1:nRows_CO,1:ncolumns_CO) !CO_JACOBIAN
            
    !jacob_all[*,nnCO[1]:(nnCO[1]+nnExt[1])-1] = ExtinctionV2_JACOBIAN 
       Jacobian_Concatenated(1:nRows_CO, &
          ncolumns_CO+1:ncolumns_CO+ncolumns_Ext) = & 
            jacobian%block(1, 7)%values(1:nRows_Ext,1:ncolumns_Ext) ! ExtinctionV2_JACOBIAN 

    !jacob_all[*,(nnCO[1]+nnExt[1]):(nnCO[1]+nnO3[1]+nnExt[1])-1] = O3_JACOBIAN
       Jacobian_Concatenated(1:nRows_CO, &
        ncolumns_CO+ncolumns_Ext+1:ncolumns_all)= & 
            jacobian%block(1, 6)%values(1:nRows_O3,1:ncolumns_O3)    ! O3_JACOBIAN 

      if (maxval(abs(Jacobian_Concatenated)) .LE. 1.E-30) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'All jacobian values are zeros. MLSCFM failed.' )
   end subroutine Copy_RadianceJacobian    !end of subroutine Copy_RadianceJacobian(...)

!
!==============
! Filter_RadianceJacobian does three things:
! (1) keep only the limb rays with mask = 0 
! (2) extract all for ptan > -2.5 (-alog10(316)).
! (3) re-arrange Jacobian
    
      
   subroutine Filter_RadianceJacobian(band9, precision9, ptanGHz,&
            Radiance_Calculated, &
            Jacobian_Concatenated , &
            Radiance_Noise, &
            Radiance_Observed, &
            Radiance_Calculated_out, &
            Jacobian_Concatenated_out , &
            Radiance_Noise_out, &
            Radiance_Observed_out)
      use VectorsModule, only: VectorValue_T
      

      ! for HDF5 files
      use MLSKINDS, only: R8, RV, RM
      
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error 
            
      use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test

      ! save to diagnosis h5 file-----------------
      
        use HDF5, only:  H5FCreate_F, H5FClose_F, HSize_T, &
              H5F_ACC_TRUNC_F, h5open_f, h5gClose_F, h5gCreate_f 
        
        use MLSHDF5, only: saveAsHDF5DS
      !-------------------------------------------
      
      type(VectorValue_T), intent(in) :: band9
      type(VectorValue_T), intent(in) :: precision9
      real(rv), dimension(:,:), pointer, intent(in) :: ptanGHz
      real(rv), dimension(:), pointer :: &
            Radiance_Calculated_out, &
            Radiance_Noise_out, &
            Radiance_Observed_out
            
      real(rv), dimension(:,:), pointer :: &
            Radiance_Calculated, &
            Jacobian_Concatenated , &
            Radiance_Noise, &
            Radiance_Observed, &
            jacob_MLS, jacob_MLS_tmp , &
            Jacobian_Concatenated_out 

      
      integer, dimension(:,:), pointer :: MASK9 => NULL()
      integer, dimension(:), pointer :: MASK9_char => NULL()
      integer, dimension(:), pointer ::   &
        ind_ptan_keep=>NULL(), mask_new=>NULL(), ind_keep=>NULL()
        

      ! Local variables
      character(len=128) :: moduleName="CFM_IO:Filter_RadianceJacobian"
      integer :: iostat, error, groupID, ierr
      real(rv), dimension(:,:), pointer :: band9_in => NULL(), &
        precision9_in => NULL()
        
        
      
      ! Local variables
      integer :: columns, ncolumns, ncolumns_mask
      integer :: rows, nrows, nrows_mask, nCh, nLimb
      integer :: i, iLimb,iCh, indNew, indOld, j
      integer :: nshape(2), nshape2(2)
      character(len=128) :: CERRMSG, chartmp
      
      

      ! save to diagnosis h5 file-----------------
      !
      !  integer :: fileID
      !
      !  
      !  call H5FCreate_F ('mlscfm_tmp.h5', H5F_ACC_TRUNC_F, fileID, iostat )
      !  if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      !  & 'Unable to open hdf5 file  mlscfm_tmph5 for output.' )
      !
      !-------------------------------------------
      
      
      ! band9 => Radiance_Observed , Radiance => Radiance_Calculated, Presision9 => Radiance_Noise
      
      ! Get band9
      
      !write (1, *) "Number of columns: ", band9%template%noChans
      !write (1, *) "Number of rows: ",    band9%template%noSurfs
      ncolumns = band9%template%noChans
      nrows = band9%template%noSurfs
      
      allocate(band9_in(nrows, ncolumns),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for band9.' )
      band9_in = transpose(reshape(band9%values(1:nrows*ncolumns, 1), &
        (/ncolumns, nrows/)))
        
      ! Get precision9
      !
      !write (1, *) "Number of columns: ", precision9%template%noChans
      !write (1, *) "Number of rows: ", precision9%template%noSurfs      
      ncolumns  = precision9%template%noChans
      nrows = precision9%template%noSurfs
      
      allocate(precision9_in(nrows, ncolumns),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for precision9.' )
      precision9_in = transpose(reshape(precision9%values(1:nrows*ncolumns, 1), &
        (/ncolumns , nrows/)))
        

      ! Code from Paul Wagner, See mail on Feb 24, 2012
      ! Add DumpMask to "use VectorsModule, only: .."  statement
      if ( .not. associated(precision9%mask) ) then
         call MLSMessage ( MLSMSG_Error, moduleName, &
         & 'Sorry--mask not associated in precision9.' )
      endif
      
      ! Get mask9
      
      !write (1, *) "Number of columns: ", precision9%template%noChans
      !write (1, *) "Number of rows: ", precision9%template%noSurfs
      
      ncolumns_mask  = precision9%template%noChans
      nrows_mask = precision9%template%noSurfs
            
      allocate(MASK9(nrows_mask, ncolumns_mask),STAT=ierr, ERRMSG=CERRMSG)
      allocate(MASK9_char(nrows_mask*ncolumns_mask),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for MASK9.' )
      DO I=1,nrows_mask*ncolumns_mask
        MASK9_char(I) = ichar(precision9%mask(I, 1))
      ENDDO
      MASK9 = transpose(reshape(MASK9_char, &
        (/ ncolumns_mask, nrows_mask/)))
      call Deallocate_test ( MASK9_char, &
            & 'Unable to deallocate array for MASK9_char', ModuleName)

!---- original IDL code ---
!; remove data for (1) mask not= 0 and (2) ptan < -2.5 
!; and stack rad etc to a vector
!; channel as the fast running indices
!
!   ind_ptan_keep = where((MASK[*,0] EQ 0) and (ptanGHz GE -2.5))
!-----------------------------
      
      nshape2 = shape(ptanGHz)
      
      if (nrows_mask .NE. nshape2(1)) then 
        call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Error! nrows_mask and nrow_ptan does not match!' )
      endif
      ind_ptan_keep => IDL_where((MASK9(1:nrows_mask,1) .EQ. 0) .and. &
        (ptanGHz(1:nrows_mask,1) .GE. -2.5))

      ! IDL_where returns [-1] if input logical array is all .FALSE.
      ! In this case we will flag "bad" and exit.
      if (ind_ptan_keep(1) .lt. 0) then
          if (maxval(ptanGHz(1:nrows_mask,1)) .LT. -2.5) then
              write(chartmp,'(E14.7)') maxval(ptanGHz(1:nrows_mask,1))
              call MLSMessage ( MLSMSG_Error, moduleName, &
                & 'maxval(ptanGHZ)='// &
                trim(chartmp)//', <-2.5, MLSCFM failed.' )
          else            
              call MLSMessage ( MLSMSG_Error, moduleName, &
                & 'Elements of MASK9(*,1) are all nonzeros, MLSCFM failed.' )
          endif
      endif
            
!---- original IDL code ---        
   ! BAND9  = reform(BAND9[ind_ptan_keep,*]) 
   !PRECISION9  = reform(PRECISION9[ind_ptan_keep,*])
   !simulatedRadiance = reform(RADIANCE[ind_ptan_keep,*])

   
!---- original IDL code ---
   !nCh = n_elements(band9[0,*])
   nshape = shape(BAND9_in) ! 125, 25
   nCh = nshape(2)
   
!---- original IDL code ---
   !nLimb = n_elements(band9[*,0])
   nLimb = size(ind_ptan_keep) 
   
!---- original IDL code ---
   !rad_MLS = fltarr(nCh*nLimb)
   allocate(Radiance_Observed_out(nCh*nLimb),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for Radiance_Observed_out.'//CERRMSG )
   Radiance_Observed_out = 0  

!---- original IDL code ---
   !noise_MLS = fltarr(nCh*nLimb)
   allocate(Radiance_Noise_out(nCh*nLimb),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for Radiance_Noise_out.'//CERRMSG )
   Radiance_Noise_out = 0
!---- original IDL code ---
   !radSim_MLS = fltarr(nCh*nLimb)
   allocate(Radiance_Calculated_out(nCh*nLimb),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for Radiance_Calculated_out.'//CERRMSG )
   Radiance_Calculated_out = 0
!---- original IDL code ---
!   i = 0
!   for iCh = 0, nCh-1 do begin
!       for iLimb = 0, nLimb-1 do begin
!           rad_MLS[i] = band9[iLimb,iCh]
!           noise_MLS[i] = precision9[iLimb,iCh]
!           radSim_MLS[i] = simulatedRadiance[iLimb,iCh]
!           i = i+1
!       endfor
!   endfor
!------------------------   
      i = 1
      do iCh = 1, nCh
       do iLimb = 1, nLimb
           Radiance_Observed_out(i) = band9_in(ind_ptan_keep(iLimb),iCh)
           Radiance_Noise_out(i) = precision9_in(ind_ptan_keep(iLimb),iCh)
           Radiance_Calculated_out(i) = &
              Radiance_Calculated(ind_ptan_keep(iLimb),iCh)
           i = i+1
       enddo
      enddo

!---- original IDL code ---
!; jacob is a bit tricky: assuming its arranged as 125 stacked [0:25] channels
!; and remove all for (1) and (2) according to mask
!   mask_new = make_array(n_elements(jacob_all[*,0]),/int, value=1)

      allocate(mask_new(size(Jacobian_Concatenated,1)),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for Radiance_Noise.'//CERRMSG )
      mask_new = 1
      i = 1

!---- original IDL code ---      
!   i = 0
!   for iLimb = 0, n_elements(mask[*,0])-1 do begin
!       for iCh = 0, n_elements(mask[0,*])-1 do begin
!           mask_new[i] = mask[iLimb,iCh]
!           if (ptanGHz[iLimb] LT -2.5) then mask_new[i] = 1 
!           i = i+1
!       endfor
!   endfor
      nshape = shape(mask9)
      do iLimb = 1, nshape(1) ! 125
        do iCh = 1, nshape(2) ! 25
           mask_new(i) = mask9(iLimb,iCh)
           if (ptanGHz(iLimb, 1) .LT. -2.5) mask_new(i) = 1 
           i = i+1
        enddo
      enddo


!---- original IDL code ---
!   ind_keep = where(mask_new EQ 0)
!   jacob_MLS_tmp = reform(jacob_all[ind_keep,*])
!; to follow Bills way, limb as the fast running index, we re-arrange Jacob
!   jacob_MLS = jacob_MLS_tmp
      nshape = shape(Jacobian_Concatenated)
      
      ind_keep => IDL_where(mask_new .EQ. 0)
      if (ind_keep(1) .lt. 0) then 
        call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'All elements in "mask_new" are nonzeros. MLSCFM failed.' )
      endif
      allocate(jacob_MLS_tmp(size(ind_keep),nshape(2)),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for jacob_MLS_tmp.'//CERRMSG )

      allocate(Jacobian_Concatenated_out(size(ind_keep),nshape(2)),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for Jacobian_Concatenated_out.'//CERRMSG )
        
      jacob_MLS_tmp = Jacobian_Concatenated(ind_keep,1:nshape(2))
!; to follow Bills way, limb as the fast running index, we re-arrange Jacob
      Jacobian_Concatenated_out = jacob_MLS_tmp

!---- original IDL code ---
!indNew = 0
!indOld = 0
!; loop over nCh chuncks in the new vector
!for i = 0, nLimb*nCh-1, nLimb do begin
!    ; loop over nLimb elements within a chunck in the new vector
!    for j = 0, nLimb-1 do begin
!        jacob_MLS[indNew,*] = jacob_MLS_tmp[indOld+j*nCh,*]
!        indNew = indNew + 1
!    endfor
!    indOld = indOld + 1
!endfor
!-----------------------------
      indNew = 1
      indOld = 1
    !; loop over nCh chuncks in the new vector
    ! nLimb = 125
    ! nCh = 25
    
      do i = 1, nLimb*nCh, nLimb 
    !    ; loop over nLimb elements within a chunck in the new vector
         do j = 0, nLimb-1
             Jacobian_Concatenated_out(indNew,1:nshape(2)) = jacob_MLS_tmp(indOld+j*nCh,1:nshape(2))
             indNew = indNew + 1
         enddo
         indOld = indOld + 1
      enddo
      ! Clean up      
      

      call Deallocate_test ( jacob_MLS_tmp, &
        & 'Unable to deallocate array for jacob_MLS_tmp', ModuleName )
      
      call Deallocate_test ( mask_new, &
        & 'Unable to deallocate array for mask_new', ModuleName )


      call Deallocate_test ( mask9, &
        & 'Unable to deallocate array for mask_new', ModuleName )

      call Deallocate_test ( band9_in, &
        & 'Unable to deallocate array for band9_in', ModuleName )

      call Deallocate_test ( precision9_in, &
        & 'Unable to deallocate array for precision9_in', ModuleName )

      call Deallocate_test ( MASK9_char, &
            & 'Unable to deallocate array for MASK9_char', ModuleName)
     
      
   end subroutine  Filter_RadianceJacobian 
   !================


   ! Write 'state;, 'radiance', 'jacobian', 'ptanGHz' to a HDF5 file, Zheng Qu
   subroutine Write_To_HDF5 (fileName, state, radiance, jacobian, newer_ptanGHz)
      ! Reference to user defined types
      use VectorsModule, only: VECTOR_T, VectorValue_T
      use MATRIXMODULE_1 , only: MATRIX_T
      ! for HDF5 files
      use MLSKINDS, only: RV
      use MLSHDF5, only: saveAsHDF5DS
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error
      use HDF5, only:  H5FClose_F, &
            H5F_ACC_TRUNC_F, H5FCreate_F

      ! Function argument type
      type(Vector_T),           intent(in) :: state
      type(Vector_T),           intent(in) :: radiance
      type(Matrix_T),           intent(in) :: jacobian
      type(VectorValue_T), intent(in) :: newer_ptanGHz
      character(len=*) , intent(in):: fileName
      
      ! Local variables
      character(len=128) :: moduleName="CFM_IO:Write_To_HDF5"
      integer :: iostat, fileID, error, groupID, ierr
      real(rv), dimension(:), pointer :: vectVals => NULL() 
      real(rv), dimension(:,:), pointer :: vectVals2 => NULL() 
      !real(rm), dimension(:,:), pointer :: matrxVals => NULL() 
      
      character(len=256)  :: a_string = 'This is a test string'
      integer, dimension(1):: Year = (/2006/) 
      
      ! Local variables
      integer :: columns, ncolumns
      integer :: rows, nrows
      character(len=128) :: CERRMSG
      
      ! initialize 
      !CALL h5open_f(error)
      !if ( error /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      !  & 'Unable to initialize hdf5 fortran interface.' )
      ! Make the Index group
      call H5FCreate_F ( trim(fileName), H5F_ACC_TRUNC_F, fileID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to open hdf5 file ' // trim(fileName) // ' for output.' )
      ! Make the Index group

      !call h5gCreate_f ( fileID, 'InputProfileData', groupID, iostat )
      !if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      !  & 'Unable to create hdf5 Index group in ' // trim(fileName) // '.' )

      

      nrows = state%quantities(1)%template%noSurfs
      allocate(vectVals(nrows),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for Temperature.'//CERRMSG )
      vectVals = state%quantities(1)%values(1:nrows, 1)
      call saveAsHDF5DS ( fileID, 'Temperature',  vectVals)
      nullify(vectVals)

      nrows = state%quantities(2)%template%noSurfs
      allocate(vectVals(nrows),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for CO2.'//CERRMSG )
      vectVals = state%quantities(2)%values(1:nrows, 1)
      call saveAsHDF5DS ( fileID, 'CO2',  vectVals)
      nullify(vectVals)
 

      nrows = state%quantities(3)%template%noSurfs
      allocate(vectVals(nrows),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for O2.'//CERRMSG )
      vectVals = state%quantities(3)%values(1:nrows, 1)
      call saveAsHDF5DS ( fileID, 'O2',  vectVals)
      nullify(vectVals)
 

 
      nrows = state%quantities(2)%template%noSurfs
      allocate(vectVals(nrows),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for CO.'//CERRMSG )
      vectVals = state%quantities(2)%values(1:nrows, 1)
      call saveAsHDF5DS ( fileID, 'CO',  vectVals)
      nullify(vectVals)

 
      nrows = state%quantities(4)%template%noSurfs
      allocate(vectVals(nrows),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for SO2.'//CERRMSG )
      vectVals = state%quantities(4)%values(1:nrows, 1)
      call saveAsHDF5DS ( fileID, 'SO2', vectVals )
      nullify(vectVals)
 

      nrows = state%quantities(5)%template%noSurfs
      allocate(vectVals(nrows),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for HNO3.'//CERRMSG )
      vectVals = state%quantities(5)%values(1:nrows, 1)
      call saveAsHDF5DS ( fileID, 'HNO3',  vectVals)
      nullify(vectVals)
 

      nrows = state%quantities(6)%template%noSurfs
      allocate(vectVals(nrows),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for O3.'//CERRMSG )
      vectVals = state%quantities(6)%values(1:nrows, 1)
      call saveAsHDF5DS ( fileID, 'O3', vectVals )
      nullify(vectVals)
      
      nrows = state%quantities(7)%template%noSurfs
      allocate(vectVals(nrows),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for ExtinctionV2R3.'//CERRMSG )
      vectVals = state%quantities(7)%values(1:nrows, 1)
      call saveAsHDF5DS ( fileID, 'ExtinctionV2R3',  vectVals)
      nullify(vectVals)
  
  

      nrows = radiance%quantities(1)%template%noChans
      ncolumns = radiance%quantities(1)%template%noSurfs

      allocate(vectVals2(nrows, ncolumns),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for Radiance.'//CERRMSG )
      vectVals2 = RESHAPE( radiance%quantities(1)%values(1:nrows*ncolumns, 1), &
        (/ nrows, ncolumns /) )
      call saveAsHDF5DS ( fileID, 'RADIANCE', vectVals2 )
      nullify(vectVals2)

      ! Save 'jacobian' values 
      nRows = jacobian%block(1, 2)%nRows
      ncolumns = jacobian%block(1, 2)%nCols
      
      allocate(vectVals2(nRows, ncolumns),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for CO_JACOBIAN.'//CERRMSG )
      vectVals2 = jacobian%block(1, 2)%values(1:nRows,1:ncolumns)
      call saveAsHDF5DS ( fileID, 'CO_JACOBIAN', vectVals2)
      nullify(vectVals2)
      
      nRows = jacobian%block(1, 6)%nRows
      ncolumns = jacobian%block(1, 6)%nCols
      
      allocate(vectVals2(nRows, ncolumns),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for O3_JACOBIAN.'//CERRMSG )
      vectVals2 = jacobian%block(1, 6)%values(1:nRows,1:ncolumns)
      call saveAsHDF5DS ( fileID, 'O3_JACOBIAN', vectVals2)
      nullify(vectVals2)
      
      nRows = jacobian%block(1, 7)%nRows
      ncolumns = jacobian%block(1, 7)%nCols
      
      allocate(vectVals2(nRows, ncolumns),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for ExtinctionV2_JACOBIAN.'//CERRMSG )
      vectVals2 = jacobian%block(1, 7)%values(1:nRows,1:ncolumns)
      call saveAsHDF5DS ( fileID, 'ExtinctionV2_JACOBIAN', vectVals2)
      nullify(vectVals2)

      
      ncolumns = newer_ptanGHz%template%noSurfs
      nrows = newer_ptanGHz%template%noChans
      allocate(vectVals2(nrows, ncolumns),STAT=ierr, ERRMSG=CERRMSG)
      if (ierr /=0) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to allocate arrays for ptanGHz.'//CERRMSG )
      vectVals2 = reshape(newer_ptanGHz%values(1:nrows*ncolumns, 1), &
        (/nrows, ncolumns/))
      
      call saveAsHDF5DS ( fileID, 'ptanGHz',vectVals2)
      nullify(vectVals2)

      call H5FClose_F ( fileID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to close hdf5 output file ' // trim(fileName) // '.' )   

        
      !call read_H5outputdata (fileName)
   end subroutine Write_To_HDF5    !end of subroutine Write_To_HDF5(...)



      subroutine Write2HDF5(outputH5fileName, &
            Radiance_Calculated_out, &
            Jacobian_Concatenated_out , &
            Radiance_Noise_out, &
            Radiance_Observed_out)
        
        ! for HDF5 files
        use MLSKINDS, only: RV
        use MLSHDF5, only: saveAsHDF5DS
        use MLSMessageModule, only: MLSMessage, MLSMSG_Error
        use HDF5, only:  H5FCreate_F, H5FClose_F, &
              H5F_ACC_TRUNC_F  
        
      real(rv), dimension(:), pointer :: &
            Radiance_Calculated_out, &
            Radiance_Noise_out, &
            Radiance_Observed_out
            
      real(rv), dimension(:,:), pointer :: &
            Jacobian_Concatenated_out 
          
        character(len=*), intent(in)  :: outputH5fileName

        ! Local variables
        character(len=128) :: moduleName="CFM_IO:Write2HDF5"
        integer :: iostat, fileID, error, ierr

        call H5FCreate_F (trim(outputH5fileName), H5F_ACC_TRUNC_F, fileID, iostat )
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to open hdf5 file ' // trim(outputH5fileName) // ' for output.' )
        call saveAsHDF5DS ( fileID, '/Radiance_Calculated',  &
            Radiance_Calculated_out)

        call saveAsHDF5DS ( fileID, '/Jacobian_Concatenated',  &
            Jacobian_Concatenated_out)

        call saveAsHDF5DS ( fileID, '/Radiance_Noise',  &
            Radiance_Noise_out)

        call saveAsHDF5DS ( fileID, '/Radiance_Observed',  &
            Radiance_Observed_out)

        call H5FClose_F ( fileID, iostat )
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Unable to close hdf5 output file ' // outputH5fileName// '.' )   

       end subroutine  Write2HDF5 !! end of subroutine Write2HDF5


       
   !================


   subroutine Read_HDF5_Spectroscopy (filename)
       use SpectroscopyCatalog_m, only: read_spectroscopy

       character(len=*), intent(in) :: filename

       call read_spectroscopy(0, filename, 'HDF5')

   end subroutine

    
    
   
    ! write data to text file--same format as pranjit did, for verification purpose
   ! Write 'state;, 'radiance', 'jacobian', 'ptanGHz' to separate text files
   subroutine Write_To_File1Txt (Temperature, CO2 , &
        O2, CO, SO2, HNO3, O3, ExtinctionV2R3, ptanGHz, &
        CO_JACOBIAN, &
        RADIANCE, O3_JACOBIAN, ExtinctionV2_JACOBIAN , &
        BAND9, PRECISION9, MASK9)
      ! Reference to user defined types
      
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error

      use MLSCommon, only: r8

      ! input variables
      
      real(r8), dimension(:), pointer ::  Temperature, CO2 , &
        O2, CO, SO2, HNO3, O3, ExtinctionV2R3
      real(r8), dimension(:,:), pointer ::  CO_JACOBIAN, &
        RADIANCE, O3_JACOBIAN, ExtinctionV2_JACOBIAN, ptanGHz, &
        BAND9, PRECISION9
      
      integer, dimension(:,:), pointer :: MASK9

      ! Local variables
      integer :: columns
      integer :: rows
      integer :: i, j
      CHARACTER(LEN=20), PARAMETER :: FMT1 = "(F16.5)"
      CHARACTER(LEN=20), PARAMETER :: FMT2 = "(E18.8)"
      character(len=*), parameter :: ModuleName = 'CFM_IO:Write_To_File1Txt'

      ! Open file to write, for 'state'
      open (unit = 1, file = "state.asc")
      !write (1, *) ""
      write (1, *) "Species 1: Temperature"
      write (1, *) "Total values: ", size(Temperature)

      do rows = 1, size(Temperature)
        write (1, FMT1, advance="no") Temperature(rows)
      end do

      !write (1, *) ""
      write (1, *) "Species 2: CO"
      write (1, *) "Total values: ", size(CO)
      do rows = 1, size(CO)
        write (1, FMT2, advance="no") CO(rows)
      end do

      !write (1, *) ""
      write (1, *) "Species 3: O2"
      write (1, *) "Total values: ", size(O2)
      do rows = 1, size(O2)
         write (1, FMT1, advance="no") O2(rows)
      end do

      !write (1, *) ""
      write (1, *) "Species 4: SO2"
      write (1, *) "Total values: ", size(SO2)
      do rows = 1, size(SO2)
         write (1, FMT2, advance="no") SO2(rows)
      end do

      !write (1, *) ""
      write (1, *) "Species 5: HNO3"
      write (1, *) "Total values: ", size(HNO3)
      do rows = 1, size(HNO3)
         write (1, FMT2, advance="no") HNO3(rows)
      end do

      !write (1, *) ""
      write (1, *) "Species 6: O3"
      write (1, *) "Total values: ", size(O3)
      do rows = 1, size(O3)
         write (1, FMT2, advance="no") O3(rows)
      end do

      !write (1, *) ""
      write (1, *) "Species 7: ExtinctionV2R3"
      write (1, *) "Total values: ", size(ExtinctionV2R3)
      do rows = 1, size(ExtinctionV2R3)
         write (1, FMT2, advance="no") ExtinctionV2R3(rows)
      end do
      ! close file opened to write
      close(1)

      ! To write 'radiance' values in a 2 dimensional format
      ! Dimension: rows = noSurfs, columns = noChans
      ! Open file to write
      open (unit = 1, file = "simulatedRadiance.asc", status='UNKNOWN')
      write (1, *) "Number of rows: ", size(Radiance,1)
      write (1, *) "Number of columns: ", size(Radiance,2)

      do rows = 1, size(radiance,2)
         do j = 1, size(Radiance,1)
             write (1, FMT1, advance="no") Radiance(j, rows)
         enddo
        write(1, *) ""
      end do
      ! close file opened to write
      close(1)

      ! Save 'jacobian' values in a file
      open (unit = 1, file = "jacobian.asc")
      write (1, *) "Species 2: CO"
      write (1, *) "Rows: ", size(CO_jacobian,1)
      write (1, *) "Columns: ", size(CO_jacobian,2)

      
      do i = 1, size(CO_jacobian,1)
            do j = 1 ,size(CO_jacobian,2)
                write (1, FMT1, advance="no") CO_jacobian(i,j)
            enddo
         write(1, *) ""
      end do

      write (1, *) "Species 6: O3"
      write (1, *) "Rows: ", size(O3_jacobian,1)
      write (1, *) "Columns: ", size(O3_jacobian,2)


      do i = 1, size(O3_jacobian,1)
            do j = 1, size(O3_jacobian,2)
               write (1, FMT1, advance="no") O3_jacobian(i,j)
            enddo
         write(1, *) ""
      end do

      write (1, *) "Species 7: ExtinctionV2"
      write (1, *) "Rows: ", size(ExtinctionV2_jacobian,1)
      write (1, *) "Columns: ", size(ExtinctionV2_jacobian,2)
      

      columns = 7 ! for EXTINCTIONV2
      do i = 1, size(ExtinctionV2_jacobian,1)
         do j = 1, size(ExtinctionV2_jacobian,2)
           write (1, FMT1, advance="no") ExtinctionV2_jacobian(i,j)
         end do
         write(1, *) ""
      end do
      ! close file opened to write
      close(1)

      ! Write ptanGHz in a text file
      ! Open file to write ptanGHz
      open (unit = 1, file = "ptanGHz.asc", status='UNKNOWN')
      write (1, *) "Number of columns: ",  size(ptanGHz,2)
      write (1, *) "Number of rows: ",     size(ptanGHz,1)
      
      do j = 1, size(ptanGHz,2)
          do i = 1, size(ptanGHz,1)
             write (1, FMT1, advance="no") ptanGHz(i,j)
          end do
      end do
      close(1)
    !============Write_to_File2===============
    
      ! Open file to write band9
      open (unit = 1, file = "band9.asc", STATUS='UNKNOWN')
      write (1, *) "Number of columns: ", size(band9,1)
      write (1, *) "Number of rows: ",    size(band9,2)

      do j = 1, size(band9,2)
          do i = 1, size(band9,1)
             write (1, FMT1, advance="no") band9(i,j)
          end do
         write(1, *) ""
      end do
      close(1)

      ! Open file to write precision9
      open (unit = 1, file = "precision9.asc", STATUS='UNKNOWN')
      write (1, *) "Number of columns: ", size(precision9,1)
      write (1, *) "Number of rows: ",    size(precision9,2)

      do j = 1, size(precision9,2)
          do i = 1, size(precision9,1)
             write (1, FMT1, advance="no") precision9(i,j)
          end do
         write(1, *) ""
      end do
      close(1)


      ! Code from Paul Wagner, See mail on Feb 24, 2012
      ! Add DumpMask to "use VectorsModule, only: .."  statement
      if ( .not. associated(mask9) ) then
         call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Sorry--mask not associated for precision9.' )   
      
      endif
      !Open file to write mask
      open (unit = 1, file = "mask9.asc", status='UNKNOWN')
      write (1, *) "Number of columns: ", size(mask9,1)
      write (1, *) "Number of rows: ", size(mask9,2)
      
      do j = 1, size(mask9,2)
          do i = 1, size(mask9,1)
             write (1, fmt='(z3.2)', advance="no") (mask9(i,j))
          end do
         write(1, *) ""
      end do
      close(1)
      
   end subroutine  Write_To_File1Txt  
   
   
   
   ! Read DACS filter shape file and add them to the DACS filter shapes
   ! database
   subroutine ReadDACSFilterShapes (fileName)
      use FilterShapes_m, only: open_filter_shapes_file, &
                                read_DACS_filter_shapes_file, &
                                close_filter_shapes_file

      character(len=*), intent(in) :: fileName
      integer :: lun, fileIndex

      ! Executables
      call open_filter_shapes_file (trim(fileName), lun, fileIndex)
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
      call open_antenna_patterns_file ( trim(fileName), lun )
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
      call open_filter_shapes_file ( trim(fileName), lun, fileIndex )
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
      call open_pointing_grid_file ( trim(fileName), lun )
      call read_pointing_grid_file ( lun, 0 )
      call close_pointing_grid_file ( lun )
   end subroutine

   ! Read PFA file and populate the PFA database
   subroutine ReadPFAFile (filename)
      use PFADataBase_m, only: process_PFA_File
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error

      character(len=*), intent(in) :: filename
      integer :: num

      num = process_PFA_file (trim(filename), 0)
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
      call mls_closeFile(file)
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
! Revision 1.11  2016/06/14 17:47:07  pwagner
! Added optional DSName arg to read+ptan
!
! Revision 1.10  2016/03/29 20:29:52  pwagner
! Restored functionality to validate_path; removed redundant ' implicit none' staements
!
! Revision 1.9  2015/08/05 20:21:29  pwagner
! Modified to compile properly with v4
!
! Revision 1.7  2011/12/15 18:27:44  honghanh
! Documentation and code clean up, including removing unused and broken
! subroutines.
!
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
