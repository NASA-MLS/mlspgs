program client
   use PVM
   use CFM
   use CFM, only: QUANTITYTEMPLATE_T
   use machine, only: getarg
   use MLSStrings, only: ReadIntsFromChars

   implicit none

!---------------------------- RCS Ident Info ------------------------------
   character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
   character (len=*), parameter :: IdParm = &
       "$Id$"
   character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

   integer :: tid, info, bufid, msgtag
   integer :: server
   character(len=10) :: input

   call PVMFMyTid(tid)

   if (tid <= 0) &
      call MLSMessage (MLSMSG_Error, moduleName, &
         "Can't contact PVM daemon")
   print *, "Client starts with tid ", tid

   call getarg(1, input)
   call ReadIntsFromChars (trim(input), server)

   print *, "server is ", server

   if (server <= 0) &
      call MLSMessage (MLSMSG_Error, moduleName, &
         "Invalid server TID")

   call PVMFInitSend(PVMDataDefault, bufid)
   if (bufid < 0) &
      call MLSMessage(MLSMSG_Error, moduleName, &
         "Can't create send buffer")

   call PVMF90Pack(1, info)
   if (info < 0) &
      call MLSMessage(MLSMSG_Error, moduleName, &
         "Can't pack signal")
   call PVMF90Pack('hello', info)
   if (info < 0) &
      call MLSMessage(MLSMSG_Error, moduleName, &
         "Can't pack signal")

   call PVMFSend(server, msgtag, info)
   if (info /= 0) &
      call MLSMessage (MLSMSG_Error, moduleName, &
         "Can't send signal to server")

   call PVMFExit(tid, info)

end program client

! $Log$
! Revision 1.1  2011/03/15 15:23:50  honghanh
! Initial imports
!
