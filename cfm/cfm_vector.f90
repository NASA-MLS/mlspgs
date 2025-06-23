! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module CFM_Vector_m
   use VectorsModule, only: Vector_T, VectorTemplate_T, VectorValue_T, &
                            DestroyVectorTemplateInfo, CloneVector, &
                            CloneVectorQuantity
   use QuantityTemplates, only: QuantityTemplate_T
   use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
         MLSMSG_Error, MLSMSG_Warning, MLSMSG_Deallocate
   use Allocate_Deallocate, only: allocate_Test, deallocate_test
   use MLSCommon, only: r8
   use String_table, only: create_string

   implicit none

   public :: CreateVector, DestroyAgileVectorContent, DestroyVectorValueContent
   public :: CreateValue4AgileVector, CreateAgileVector, AddValue2Vector
   public :: CloneAgileVector

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

    private

    contains

    ! This method is obsolete.
    type(Vector_T) function CreateVector (vectorTemplate, qtyDatabase, name) &
             result (vector)
        ! template listing the quantities to be stored in this vector
        type (VectorTemplate_T), intent (in), target :: vectorTemplate
        ! quantity template database to retrieve the template for the quantities
        ! to be stored in this vector
        type (QuantityTemplate_T), dimension(:), intent(in), target :: qtyDatabase
        ! name of the vector as string
        character(len=*), optional :: name

        integer :: quantity
        integer :: status

        vector%template = vectorTemplate
        allocate (vector%quantities(vectorTemplate%noQuantities), stat=status)
        if (status /= 0) call MLSMessage (MLSMSG_Error, moduleName, &
             MLSMSG_Allocate // "Vector quantities")

        do quantity = 1, vectorTemplate%noQuantities
            vector%quantities(quantity)%index = quantity
            vector%quantities(quantity)%template = &
                qtyDatabase(vectorTemplate%quantities(quantity))
        end do
        call CreateValues(vector)

        if (present(name)) then
            vector%name = create_string(name)
        else
            vector%name = 0
        end if
    end function

    ! Initializes a vector to which vector values can be dynamically added.
    type(Vector_T) function CreateAgileVector (name) result (vector)
        character(len=*), optional :: name

        ! this vector's template has no quantities array
        vector%template%name = 0
        vector%template%noQuantities = 0
        vector%template%totalInstances = 0
        vector%template%totalElements = 0
        nullify(vector%template%quantities)
        nullify(vector%quantities)

        if (present(name)) then
            vector%name = create_string(name)
        else
            vector%name = 0
        end if
    end function

    ! Clone a vector content, including its values
    ! but the values' templates are only cloned as references level.
    subroutine CloneAgileVector (y, x)
        type(Vector_T), intent(in) :: x
        type(Vector_T), intent(out) :: y

        integer :: i

        call CloneVector(y, x)
        do i = 1, size(y%quantities)
            call CloneVectorQuantity(y%quantities(i), x%quantities(i))
        end do
    end subroutine

    ! Cleans up objects that are allocated inside an agile vector.
    subroutine DestroyAgileVectorContent (v)
        type(Vector_T), intent(inout) :: v
        integer :: i

        do i = 1, size(v%quantities)
            call DestroyVectorValueContent(v%quantities(i))
        end do
        call DestroyVectorTemplateInfo(v%template)
    end subroutine

    ! Only use this if the vector value does not belong to
    ! any vector.
    ! To destroy an entire vector, use the methods dedicated
    ! to destroy vector for better efficiency.
    subroutine DestroyVectorValueContent (vv)
        type(VectorValue_T), intent(inout) :: vv

        call deallocate_test (vv%values, "vv%values", moduleName)
        call deallocate_test ( vv%mask, "vv%mask", ModuleName )
    end subroutine

    ! Add a VectorValue_T object to a Vector_T object.
    ! The vector has to be an agile vector.
    ! vectorvalue doesn't have to be unique.
    ! Value will be added to the end of vector.
    ! Note that memory allocated to vectorvalue will be
    ! reused by this vector, so don't add the same
    ! vectorvalue to 2 different vectors.
    subroutine AddValue2Vector (vector, vectorvalue)
        type(Vector_T), intent(inout) :: vector
        type(VectorValue_T), intent(in) :: vectorvalue

        type(VectorValue_T), dimension(:), pointer :: temp
        integer :: status

        vector%template%noQuantities = vector%template%noquantities + 1
        temp => vector%quantities

        allocate(vector%quantities(vector%template%noquantities), stat=status)
        if (status /= 0) call MLSMessage (MLSMSG_Error, moduleName, &
            MLSMSG_Allocate // "AddValue2Vector")
        if (associated(temp)) vector%quantities(1:(vector%template%noquantities-1)) = temp

        vector%quantities(vector%template%noquantities) = vectorvalue
        vector%quantities(vector%template%noquantities)%index = vector%template%noquantities

        vector%template%totalInstances = vector%template%totalinstances + &
        vectorvalue%template%noinstances
        vector%template%totalElements = vector%template%totalelements + &
                                        vectorvalue%template%noinstances * &
                                        vectorvalue%template%instanceLen
    end subroutine

    ! Creates a VectorValue_T object to add to an agile vector.
    ! The value can be empty, which means it only contains a QuantityTemplate_T.
    ! Or it can be filled with explicit values given by the value array.
    ! Or it can be filled with the same value given by a spreadvalue number.
    ! Value can also be accompanied by a mask to distinguish good values 
    ! and bad values.
    type(VectorValue_T) function CreateValue4AgileVector (template, value, spreadvalue, &
    mask) result(vectorvalue)
        use MLSStrings, only: writeIntsToChars
        use VectorsModule, only: remapVectorValue, remapVectorMask

        type(QuantityTemplate_T), intent(in) :: template
        real(r8), dimension(:), intent(in), optional :: value
        real(r8), intent(in), optional :: spreadvalue
        character, dimension(:), intent(in), optional :: mask

        integer :: row, col
        character(len=10) :: int1 = "          ", int2 = "          "
        integer :: valueSize

        row = template%nochans * template%nosurfs
        col = template%noinstances

        vectorvalue%template = template

        valueSize = max( 1, vectorvalue%template%noChans * &
           & vectorvalue%template%noSurfs * &
           & vectorvalue%template%noInstances * &
           & vectorvalue%template%noCrossTrack )
        if (.not. present(value) .and. .not. present(spreadvalue)) then
            call allocate_test(vectorvalue%value1, valueSize, "vectorvalue%values", modulename)
            vectorvalue%value1 = 0.0_r8
            call remapVectorValue ( vectorvalue )
            return
        end if

        if (present(value) .and. present(spreadvalue)) &
        call MLSMessage(MLSMSG_Error, modulename, "Cannot specify both value and spreadvalue.")

        ! Fortran do not have lazy evaluation of conditional clause
        if (present(value)) then
            if (row * col /= size(value)) then
                call writeintstochars(row * col, int1)
                call writeintstochars(size(value), int2)

                call MLSMessage (MLSMSG_Error, moduleName, &
                "Incorrect size, expect " // int1 // " elements, got " // int2)
            end if
        end if

        ! This simply doesn't work any longer after we moved to remapping
        if ( .false. ) then
          call allocate_test(vectorvalue%values, row, col, "vectorvalue%values", modulename)
          if (present(spreadvalue)) then
              vectorvalue%values = spreadvalue
          else
              vectorvalue%values = reshape(value, shape(vectorvalue%values))
          end if
        else
          call allocate_test ( vectorvalue%value1, &
            & valueSize, &
            & trim(vectorvalue%AllocationName), moduleName )
          if (present(spreadvalue)) then
              vectorvalue%value1 = spreadvalue
          else
              vectorvalue%value1 = value
          end if

          call remapVectorValue ( vectorvalue )
        endif

        if (present(mask)) then
            if (size(vectorvalue%values) /= size(mask)) &
            call MLSMessage (MLSMSG_Error, modulename, "Size of mask differs from shape of value.")
            call allocate_test(vectorvalue%mask1, valueSize, "vectorvalue%mask", modulename)
            vectorvalue%mask1 = mask
          call remapVectorMask ( vectorvalue )
        end if
    end function CreateValue4AgileVector

    ! =====     Private Procedures     =====================================
    subroutine CreateValues ( Vector)
        ! Allocate space for the values of a vector.
        type(Vector_T), intent(inout) :: Vector
        integer :: QTY
        real(r8), parameter :: MYHUGE = 1.0e15

        do qty = 1, size(vector%quantities)
            call allocate_test ( vector%quantities(qty)%values, &
                & vector%quantities(qty)%template%noChans * &
                & vector%quantities(qty)%template%noSurfs, &
                & vector%quantities(qty)%template%noInstances, &
                & "vector%quantities(qty)%%values", ModuleName )
            vector%quantities(qty)%values = 0.0_r8
        end do
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
! Revision 1.14  2011/12/25 03:39:06  honghanh
! Fix cloning and memory cleanup method
! in forwardmodel2 subroutines.
!
! Revision 1.13  2011/12/15 18:27:45  honghanh
! Documentation and code clean up, including removing unused and broken
! subroutines.
!
! Revision 1.12  2011/11/02 03:50:33  honghanh
! Add code to clean up values array in DestroyVectorValueContent
!
! Revision 1.11  2011/11/01 22:16:11  honghanh
! Add API to destroy individual QuantityTemplate_T and VectorValue_T
!
! Revision 1.10  2011/10/31 19:12:32  honghanh
! Add "intent(in)" to input in CreateValue4AgileVector
!
! Revision 1.9  2011/10/19 19:33:10  honghanh
! Adding DestroyAgileVectorContent to CFM_Vector_m
!
! Revision 1.8  2011/03/24 15:16:46  honghanh
! Add new interfaces for creating vector and vector values without going through quantity template databases
!
! Revision 1.7  2010/11/03 20:17:01  honghanh
! Add name as an optional argument to CreateVector.
!
! Revision 1.6  2010/06/29 16:40:23  honghanh
! Remove all function/subroutine and user type forwarding from
! all CFM modules except for from cfm.f90
!
! Revision 1.5  2010/06/29 15:53:45  honghanh
! Add copyright comments and support for CVS log in the file
!
