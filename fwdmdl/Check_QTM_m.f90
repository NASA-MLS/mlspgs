! Copyright 2015, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Check_QTM_m
!=============================================================================


  implicit NONE

  private

  public :: Check_QTM

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! --------------------------------------------------  Check_QTM  -----
  subroutine Check_QTM ( Grids_Tmp, Grids_f, QTM_HGrid, UsingQTM )
    ! Check whether temperature and all species either all have QTM HGrid,
    ! or none do.  If they all do, set UsingQTM and find the QTM HGrid with
    ! the finest resolution and associate it with QTM_HGrid.

    use HGridsDatabase, only: HGrid_T
    use Intrinsic, only: Lit_Indices, L_Temperature
    use Load_SPS_Data_M, only: Grids_T
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Output_m, only: Output
    use String_Table, only: Display_String
    use Toggles, only: Emit, Levels, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use VectorsModule, only: VectorValue_T

    type (Grids_T), intent(in) :: Grids_Tmp  ! All the coordinates for Temperature
    type (Grids_T), intent(in) :: Grids_f    ! All the coordinates for VMR
    type (HGrid_T), pointer, intent(out) :: QTM_HGrid ! HGrid that has finest QTM resolution.
    logical, intent(out) :: UsingQTM         ! Temperature and all species have QTM hGrids

    integer :: I, K     ! Loop indices
    integer :: Me = -1  ! String index for trace
    integer :: Pos
    logical :: QTM_fail
    logical :: SpsQTM   ! Species being examined in Grids_F has QTM hGrid
    integer, parameter :: Width = 100 ! of line for list of quantities

    call trace_begin ( me, 'ForwardModel.Check_QTM', &
      & cond=toggle(emit) .and. levels(emit) > 0  )

    QTM_fail = .false.  ! Assume no errors

    ! Verify that temperature and all the species either all have QTM hGrids
    ! or none do.  If they all do, set QTM_HGrid to the one with the finest
    ! resolution.
    QTM_hGrid => grids_tmp%qtyStuff(1)%qty%template%the_hGrid
    usingQTM = grids_tmp%isQTM(1)
    do k = 1, size(grids_f%qtyStuff)
      spsQTM = grids_f%isQTM(k)
      QTM_fail = QTM_fail .or. ( spsQTM .neqv. usingQTM )
      if ( .not. QTM_fail ) then
        if ( QTM_hGrid%QTM_tree%level < &
           & grids_f%qtyStuff(k)%qty%template%the_hGrid%QTM_tree%level ) &
           & QTM_hGrid => grids_f%qtyStuff(k)%qty%template%the_hGrid
      end if
    end do
    if ( QTM_fail ) then
      do i = 1, 2 ! 1 = with QTM, 2 = without QTM
        call output ( 'Quantities with' // trim(merge('   ','out',i==1)) // &
                    & ' QTM HGrids:', advance='yes' )
        pos = 0
        if ( i == 1 .eqv. usingQTM ) call display_string ( &
          & lit_indices(l_temperature), before=' ', pos=pos, width=width )
        do k = 1, size(grids_f%qtyStuff)
          spsQTM = grids_f%isQTM(k)
          if ( i == 1 .eqv. spsQTM ) call display_string ( &
            & lit_indices(grids_f%qtyStuff(k)%qty%template%quantityType), &
            & before=' ', pos=pos, width=width )
        end do
        call newLine
      end do
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Some quantities have QTM hGrids but some do not.' )
    end if

    if ( .not. usingQTM ) nullify ( QTM_HGrid )

    call trace_end ( 'ForwardModel.Both_Sidebands_Setup', &
      & cond=toggle(emit) .and. levels(emit) > 0  )

  end subroutine Check_QTM

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Check_QTM_m

! $Log$
! Revision 2.1  2016/09/12 23:49:51  vsnyder
! Initial commit
!
