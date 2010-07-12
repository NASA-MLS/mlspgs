! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
program test_addtodatabase
   use CFM, only: QuantityTemplate_T, MLSFile_T, &
                  AddQuantityTemplateToDatabase, &
                  AddFileToDatabase

   implicit none

!---------------------------- RCS Ident Info ------------------------------
   character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
   character (len=*), parameter :: IdParm = &
       "$Id$"
   character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

   type(QuantityTemplate_T), dimension(:), pointer :: quantities
   type(MLSFile_T), dimension(:), pointer :: files
   type(MLSFile_T), target :: f1, f2
   type(QuantityTemplate_T), target :: q1, q2
   integer :: i, oldsize, newsize
   character(len=3) :: answer

   print *, ""
   print *, "-----------------"
   print *, "SUCCESS CRITERIA:"
   print *, "-----------------"
   print *, "1. Size of database increase by 1 after adding an item."
   print *, "2. The element at the returned index into the database is the added item."
   print *, "3. The returned index must be greater than 0 but less than or equal to "
   print *, "   the size of the database (or the returned index is valid)."
   print *, "4. Can add to an empty database."
   print *, "5. Can add to database with existing items."
   print *, "6. Can add duplicate items."

   print *, ""
   print *, "*******************************************************************************"
   print *, ""

   nullify(quantities, files)

   print *, "-----------------------------------"
   print *, "Test AddQuantityTemplateToDatabase:"
   print *, "-----------------------------------"
   q1%name = 1
   q2%name = 2
   print *, "1. Add to an empty database:"
   i = AddQuantityTemplateToDatabase(quantities, q1)
   if (size(quantities) == 1) then
      answer="YES"
   else
      answer="NO "
   end if
   print *, "- Is size of database 1? ", answer
   if (i == 1) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is returned index 1? ", answer
   if (quantities(1)%name == 1) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is element at index 1 the added item? ", answer
   print *, "2. Add to a database with existing item(s):"
   oldsize = size(quantities)
   i = AddQuantityTemplateToDatabase(quantities, q2)
   newsize = size(quantities)
   if (newsize == (oldsize + 1)) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Does size increase by 1? ", answer
   if (i <= newsize .and. i > 0) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is the returned index valid? ", answer
   if (quantities(i)%name == 2) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is element at the returned index the added item? ", answer
   print *, "3. Add duplicate item:"
   oldsize = size(quantities)
   i = AddQuantityTemplateToDatabase(quantities, q2)
   newsize = size(quantities)
   if (newsize == (oldsize + 1)) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Does size increase by 1? ", answer
   if (i <= newsize .and. i > 0) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is the returned index valid? ", answer
   if (quantities(i)%name == 2) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is element at the returned index the added item? ", answer
   i = 0
   do newsize = 1, size(quantities)
      if (quantities(newsize)%name == 2) i = i + 1
   end do
   if (i == 2) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is the number of duplicated items 2? ", answer

   print *, ""
   print *, "*******************************************************************************"
   print *, ""

   print *, "-----------------------"
   print *, "Test AddFileToDatabase:"
   print *, "-----------------------"
   f1%name = "file1"
   f2%name = "file2"
   print *, "1. Add to an empty database:"
   i = AddFileToDatabase(files, f1)
   if (size(files) == 1) then
      answer="YES"
   else
      answer="NO "
   end if
   print *, "- Is size of database 1? ", answer
   if (i == 1) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is returned index 1? ", answer
   if (files(1)%name == "file1") then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is element at index 1 the added item? ", answer
   print *, "2. Add to a database with existing item(s):"
   oldsize = size(files)
   i = AddFileToDatabase(files, f2)
   newsize = size(files)
   if (newsize == (oldsize + 1)) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Does size increase by 1? ", answer
   if (i <= newsize .and. i > 0) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is the returned index valid? ", answer
   if (files(i)%name == "file2") then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is element at the returned index the added item? ", answer
   print *, "3. Add duplicate item:"
   oldsize = size(files)
   i = AddFileToDatabase(files, f2)
   newsize = size(files)
   if (newsize == (oldsize + 1)) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Does size increase by 1? ", answer
   if (i <= newsize .and. i > 0) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is the returned index valid? ", answer
   if (files(i)%name == "file2") then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is element at the returned index the added item? ", answer
   i = 0
   do newsize = 1, size(files)
      if (files(newsize)%name == "file2") i = i + 1
   end do
   if (i == 2) then
      answer = "YES"
   else
      answer = "NO "
   end if
   print *, "- Is the number of duplicated items 2? ", answer

   print *, ""
   print *, "*******************************************************************************"
   print *, ""

   print *, "----------------------------------------------------------------------------"
   print *, "VERIFICATIONS: The test program ends normally and all the answers are 'YES'."
   print *, "----------------------------------------------------------------------------"
   print *, ""
end program test_addtodatabase

! $Log$
