module PATH_ENTITIES_M
  use MLSCommon, only: I4, R4, R8
  implicit NONE
  public
! This file contains the definistion of all the PATH related entities
! Written  by: Z. Shippony, Aug/28/2000
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName = &
       "$RCSfile$"
!---------------------------------------------------------------------------

 type PATH_INDEX
   Integer(i4) :: break_point_index
   Integer(i4) :: total_number_of_elements
 end type PATH_INDEX

 type PATH_VECTOR
   Real(r8), DIMENSION(:), POINTER :: values => NULL()
 end type PATH_VECTOR

 type PATH_VECTOR_2D
   Real(r8), DIMENSION(:,:), POINTER :: values => NULL()
 end type PATH_VECTOR_2D

 type PATH_DERIVATIVE
   Real(r4), DIMENSION(:,:,:), POINTER :: values => NULL()   ! (various,mnp,mxco)
 end type PATH_DERIVATIVE

 type PATH_BETA
   Real(r8), DIMENSION(:), POINTER :: values => NULL()
   Real(r8), DIMENSION(:), POINTER :: t_power => NULL()
   Real(r8), DIMENSION(:), POINTER :: dbeta_dw => NULL()
   Real(r8), DIMENSION(:), POINTER :: dbeta_dn => NULL()
   Real(r8), DIMENSION(:), POINTER :: dbeta_dnu => NULL()
 end type PATH_BETA

!type PATH_CROSS_SECTION
!  Real(r8), DIMENSION(:), POINTER :: values
!  Real(r8), DIMENSION(:), POINTER :: t_power
!end type PATH_CROSS_SECTION

end module PATH_ENTITIES_M
! $Log$
! Revision 1.4  2001/03/26 17:56:14  zvi
! New codes to deal with dh_dt_path issue.. now being computed on the fly
!
! Revision 1.3  2001/03/09 23:57:02  vsnyder
! Declare pointer components to be initially nullified.
!
! Revision 1.2  2001/01/31 01:08:48  zvi
! New version of forward model
!
! First version 2000/08/28 21:56:12  zvi
