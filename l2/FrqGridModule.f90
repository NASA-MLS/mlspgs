! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module FrqGridModule
  !=============================================================================

  implicit none

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  public :: FrqGrid

contains

  subroutine FrqGrid ( Root, VectorDatabase, configDatabase )

  ! THIS ROUTINE IS INCOMPLETE.  IT IS MEANT TO BE A PLACE FOR BILL TO
  ! HANG A REFERENCE TO FREQUENCY-GRID CALCULATION SOFTWARE.  EVENTUALLY,
  ! WE WANT THAT SOFTWARE TO USE THE SAME INFRASTRUCTURE AS THE REST OF
  ! THE SOFTWARE.  AT PRESENT, THIS ROUTINE ISN'T CALLED.  THE IDEA IS
  ! TO CALL IT AFTER THE NECESSARY PARTS OF THE L2CF HAVE BEEN PROCESSED.
  ! MAYBE THE RETRIEVE SECTION IS THE RIGHT PLACE.  ANYWAY, AT PRESENT,
  ! INIT_TABLES_MODULE DOESN'T ALLOW THE FREQUENCYGRID SPECIFICATION TO
  ! APPEAR ANYWHERE.

  ! Process the FrequencyGrid specification.  It has an 'Atmos' field
  ! that is required to be a vector, and a 'Frequencies' field that is
  ! required to be a numeric array.  The values of the Frequencies field
  ! must have either no units, or frequency units, and at least one of them
  ! has to have frequency units.

  ! "Processing" the FrequencyGrid specification results in computing a
  ! frequency grid.

    use Allocate_Deallocate, only: Allocate_Test
    use Expr_M, only: Expr
    use ForwardModelConfig, only: ForwardModelConfig_T
    use Init_Tables_Module, only: F_Atmos, F_Frequencies
    use Intrinsic, only: PHYQ_Dimensionless, PHYQ_Frequency
    use MLSCommon, only: RK => R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MoreTree, only: Get_Field_Id
    use Tree, only: Decoration, Nsons, Source_Ref, Subtree
    use VectorsModule, only: Vector_T

    ! Dummy arguments:
    integer, intent(in) :: Root         ! Of the relevant subtree of the AST
    ! Indexes an n_cf vertex
    type(vector_T), dimension(:), intent(inout), target :: VectorDatabase
    type(forwardModelConfig_T), dimension(:), pointer :: configDatabase

    ! Local variables:
    type(vector_T), pointer :: Atmos    ! Atmospheric state vector
    integer :: Error
    integer :: Field                    ! Of the FrqGrid specification
    integer :: I, J                     ! Subscript, loop inductor
    integer :: Son                      ! of Root
    real(rk), pointer :: TheFrequencies(:) ! from the Frequencies field
    integer :: TheUnits                 ! of the Frequencies field
    integer :: Type                     ! Type of value returned by EXPR
    integer :: Units(2)                 ! Units of value returned by EXPR
    double precision :: Value(2)        ! Value returned by EXPR

    ! Error codes
    integer, parameter :: NoUnits = 1   ! No field has frequency units
    integer, parameter :: WrongUnits = NoUnits + 1 ! not freq and not dimless

    error = 0
    do i = 2, nsons(root)
      son = subtree(i,root)
      field = get_field_id(son)
      select case ( field )
      case ( f_atmos )
        atmos => vectorDatabase(decoration(decoration(subtree(2,son))))
      case ( f_frequencies )
        theUnits = phyq_dimensionless
        call allocate_test ( theFrequencies, nsons(root)-1, 'TheFrequencies', &
          & moduleName )
        do j = 2, nsons(root)
          call expr ( subtree(j,son), units, value, type )
          ! Make sure the value has frequency units.  Dimensionless is OK,
          ! so long as at least one has frequency units.
          if ( units(1) == phyq_frequency ) then
            theUnits = phyq_frequency
          else if ( units(1) /= phyq_dimensionless ) then
            call announceError ( wrongUnits, subtree(j,son) )
          end if
          theFrequencies(j-1) = value(1)
        end do
        if ( theUnits /= phyq_frequency ) then
          call announceError ( noUnits, subtree(j,son) )
        end if
      end select
    end do ! i = 2, nsons(root)

    call MLSMessage ( MLSMSG_Error, moduleName, &
      & "No frequency grid computed -- error in configuration" )

    ! Bill: Put the call(s) to the frequency-gridding procedures here.
    ! The atmospheric state vector is "Atmos".  The frequencies are
    ! "TheFrequencies".  They've been checked to make sure they have
    ! frequency units.  It's OK in the input if some have no units,
    ! but at least one has to have frequency units.  They're scaled
    ! to the "standard" frequency unit, which is MHz.  Let me know
    ! if you need extra stuff, say an output vector.

  contains

    subroutine AnnounceError ( Why, Where )

      use Lexer_Core, only: Print_Source
      use Output_M, only: Output

      integer, intent(in) :: Why
      integer, intent(in) :: Where

      error = max(error,1)
      call output ( '***** At ' )
      call print_source ( source_ref(where) )
      call output ( ', RetrievalModule complained: ' )
      select case ( why )
      case ( noUnits )
        call output ( 'No field has frequency units', advance='yes' )
      case ( wrongUnits )
        call output ( 'Field value must be unitless or have frequency units', &
          & advance='yes' )
      end select
    end subroutine AnnounceError
  end subroutine FrqGrid

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module FrqGridModule

! $Log$
! Revision 2.1  2003/04/30 00:09:25  vsnyder
! Initial commit
!
