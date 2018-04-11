! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Spectroscopy_Types

  use Intrinsic, only: L_None
  use MLSKINDS, only: R8
  implicit NONE
  private

  ! Public procedure
  public :: AddLineToDatabase

  ! Parameters
  integer, public, parameter :: MaxContinuum = 6

  ! Public types:
  type, public :: Line_T           ! One line in the spectrum for a species
    integer :: Line_Name = 0       ! Sub_rosa index
    real(r8) :: DELTA              ! Delta interference coefficient at 300K 1/mb
    real(r8) :: EL                 ! Lower state energy cm-1
    real(r8) :: GAMMA              ! Gamma interference coefficient at 300K 1/mb
    real(r8) :: N                  ! Temperature power dependence of w
    real(r8) :: N1                 ! Temperature dependency of delta
    real(r8) :: N2                 ! Temperature dependency of gamma
    real(r8) :: NS                 ! Pressure shift dependency on temperature
    real(r8) :: PS                 ! Pressure shift parameter in MHz/mbar
    real(r8) :: STR                ! Integrated spectral intensity
                                   ! Log(nm**2 MHz) at 300 K
    real(r8) :: V0                 ! Line center frequency MHz
    real(r8) :: W                  ! Collision broadening parameter
                                   ! MHz/mbar at 300 K
    logical :: UseYi               ! Delta /= 0 .or. Gamma /= 0
    integer, dimension(:), pointer :: QN=>NULL()  ! Optional quantum numbers.
                                   ! The first element is assumed to be the
                                   ! quantum number format;  The low-order digit
                                   ! of the quantum number format is the number
                                   ! of quantum number pairs.  The rest,
                                   ! necessarily an even number, are the actual
                                   ! quantum numbers.  The extent of this array
                                   ! is twice the low-order digit of the format,
                                   ! plus one.
    integer, dimension(:), pointer :: Signals=>NULL()   ! List of signal indices for line
    integer, dimension(:), pointer :: Sidebands=>NULL() ! Sidebands for above bands (-1,0,1)
    logical, dimension(:), pointer :: Polarized=>NULL() ! Process this signal and
                                   ! sideband using the polarized model
  end type Line_T

  type, public :: Catalog_T        ! Catalog entry for a species
    real(r8) :: Continuum(MaxContinuum)  ! Continuum coefficients
    real :: DefaultIsotopeRatio = 1.0
    integer, pointer :: Lines(:)=>NULL() ! Indices in Lines database
    real(r8) :: Mass               ! Molecular mass in AMU
    integer :: Molecule = l_none   ! L_...; l_none => no catalog entry
    logical, pointer :: Polarized(:)=>NULL() ! Used only in catalog extract to
                                   ! indicate that the lines(:) are to be
                                   ! processed with the polarized model
    real(r8) :: Qlog(3)            ! Logarithm of the partition function
                                   ! at 300, 225, and 150 K
    integer :: Species_Name = 0    ! Sub_rosa index of s_spectra label
  end type Catalog_T

  ! The lines database:
  type(Line_T), pointer, public, save :: Lines(:) => NULL()

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ------------------------------------------  AddLineToDatabase  -----
  integer function AddLineToDatabase ( Database, Item )
  ! Add a line to the Lines database, creating the database
  ! if necessary.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type(line_T), pointer, dimension(:) :: Database
    type(line_T), intent(in) :: Item


    ! Local variables
    type(line_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddLineToDatabase = newSize
  end function AddLineToDatabase

! ------------------------------------------------  not_used_here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Spectroscopy_Types

! $Log$
! Revision 2.7  2018/04/11 22:25:23  vsnyder
! Remove USE for unused names
!
! Revision 2.6  2015/03/28 02:06:20  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.5  2014/09/05 20:54:24  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.4  2014/04/04 19:22:33  vsnyder
! More comments about QN component of Lines_t
!
! Revision 2.3  2011/11/08 19:49:25  vsnyder
! Get L_None from Intrinsic module
!
! Revision 2.2  2011/11/01 22:15:41  vsnyder
! Add a comment to describe the Lines pointer
!
! Revision 2.1  2011/11/01 20:50:06  vsnyder
! Initial commit
!
