! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

Module Bill_GasAbsorption

   USE MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error, &
     & MLSMSG_DeAllocate
  implicit NONE
  private
  public :: get_beta_bill

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
  "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

  SUBROUTINE get_beta_bill (T, PB, F, RH, VMR_in, NS, ABSC, Catalog, LosVel )

!==============================================================
!      CALCULATE CLEAR-SKY ABSORPTION COEFFICIENT AT F AND T
!      USING BILL'S NEW SPECTROSCOPY DATA. 
!      LATEST UPDATE: J.JIANG, NOVEMBER 14, 2001
!==============================================================

    use Get_Beta_Path_m, only: CREATE_BETA
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT, ALLOCATEONESLABS, &
      & DESTROYCOMPLETESLABS
    use MLSCommon, only: R8, RP, IP
    use Molecules, only: L_H2O, L_H2O_18, L_HNO3, L_N2, L_N2O, L_O_18_O, &
      & L_O2, L_O3
    use Physics, only: H_OVER_K
    use SLABS_SW_M, only: GET_GL_SLABS_ARRAYS
    use SpectroscopyCatalog_m, only: CATALOG_T, LINES
    use WaterVapor, only: RHtoEV

    !-----------------
    ! INPUTS
    !-----------------

    INTEGER, INTENT(IN) :: NS
    REAL(r8), INTENT(IN) :: F               ! FREQUENCY IN GHz
    REAL(r8), INTENT(IN) :: T               ! TEMPERATURE (K)
    REAL(r8), INTENT(IN) :: PB              ! TOTAL AIR PRESSURE (hPa)
    REAL(r8), INTENT(IN) :: VMR_in(NS-1)    ! VMR OF SPECIES 
    REAL(r8), INTENT(IN) :: RH              ! H2O VOLUME MIXING RATIO
    REAL(rp), INTENT(IN) :: LosVel          !     OR RELATIVE HUMIDITY
    
    Type(Catalog_T), INTENT(IN) :: Catalog(:)
    Type (slabs_struct), DIMENSION(:,:), POINTER :: gl_slabs

    !------------------
    ! OUTPUTS
    !------------------
    REAL(r8), INTENT(out) :: ABSC           ! ABSORPTION COEFFICIENT (1/m)

    !------------------
    ! LOCAL VARIABLES
    !-----------------

    REAL(r8) :: B                           ! BETA (1/m/ppv)                 
    REAL(r8) :: FF                          ! FREQUENCY IN MHz               
    REAL(rp) :: myLosVel
    REAL(r8) :: P                           ! DRY AIR PARTIAL PRESSURE (hPa) 
    REAL(r8) :: VMR                         ! VOLUME MIXING RATIO            
    REAL(r8) :: VMR_H2O                     ! H2O VOLUME MIXING RATIO        
    REAL(r8) :: VMR_O3                      ! H2O VOLUME MIXING RATIO        
    REAL(r8) :: VMR_O2                      ! O2    VOLUME MIXING RATIO      
    REAL(r8) :: VMR_O_18_O                  ! O18O  VOLUME MIXING RATIO      
    REAL(r8) :: VMR_H2O_18                  ! H2O18 VOLUME MIXING RATIO      
    REAL(r8) :: VMR_N2                      ! N2 VOLUME MIXING RATIO         
    REAL(r8) :: VMR_N2O                     ! N2O VOLUME MIXING RATIO        
    REAL(r8) :: VMR_HNO3                    ! HNO3 VOLUME MIXING RATIO       
    REAL(r8) :: VP                          ! VAPOR PARTIAL PRESSURE (hPa)

    Integer(ip) :: n_sps, i, j, no_of_lines, n_ele
    integer(ip), parameter :: IPSD=1000, NU=16, NUA=8, NAB=50, NR=40, NC=2
    Integer(ip) :: status
    REAL(rp) :: bb, del_temp, cld_ext, WC(2), tanh1
    real(r8), parameter :: Boltzmhz = 1 / h_over_k
    real(rp), allocatable, dimension(:) :: PP, TT
    logical :: Do_1D, Incl_Cld
    LOGICAL, ALLOCATABLE, dimension(:) :: true_path_flags

!-----------------------------------------------------------------------------
    WC= 0._r8
    IF (RH /= 100.0_r8) THEN
       VMR_H2O = RH                     ! PH HERE IS WATER VAPOR MIXING RATIO
       VP=VMR_H2O*PB                    ! VP IS VAPOR PRESSURE, PB IS TOTAL
       P=PB-VP                          ! PRESSURE, P IS DRY-AIR PRESSURE
    ELSE IF (VMR_H2O == 100.0_r8) THEN
       CALL RHtoEV(T, 100.0_r8, VP)     ! RH HERE IS 100% RELATIVE HUMIDITY 
       P = PB-VP
       VMR_H2O = VP/(max(1.e-19_r8, P))
    END IF

    VMR_O2     = 0.2095_r8
    VMR_N2     = 0.805_r8
    VMR_O_18_O = VMR_O2*0.00409524_r8 
    VMR_H2O_18 = VMR_H2O*0.00204_r8
    VMR_O3     = VMR_in(1)
    VMR_N2O    = VMR_in(2)
    VMR_HNO3   = VMR_in(3)

    B = 0.0_r8
    FF = F*1000.0_r8

    n_sps = Size(Catalog)
!    maxVert = Ng +1
!    n_ele = 2*maxVert
    n_ele = 1  ! number of pressure levels along the path, in our case =1

    allocate ( true_path_flags(n_ele), stat=status )
    true_path_flags = .true.
    allocate ( gl_slabs(n_ele,n_sps), stat=status )
    if ( status /= 0 ) &
      & CALL MLSMessage(MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // ' gl_slabs ')
    allocate ( pp(n_ele), stat=status )
    if ( status /= 0 ) &
      & CALL MLSMessage(MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // ' pp ')
    allocate ( tt(n_ele), stat=status )
    if ( status /= 0 ) &
      & CALL MLSMessage(MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // ' tt ')
    do i = 1, n_sps
      no_of_lines =  size(Catalog(i)%Lines)
      do j = 1, n_ele
          Call AllocateOneSlabs ( gl_slabs(j, i), no_of_lines )
      enddo
    enddo

    pp(1) = p
    tt(1) = t
                              ! Bill uses km/sec 
    myLosVel=losVel*0.00_rp   ! The Doppler correction already been done 
                              ! in the FullCloudForwardModel, so set it 0

    call get_gl_slabs_arrays ( Catalog, PP(1:n_ele), TT(1:n_ele), myLosVel, &
      & gl_slabs, Do_1D, true_path_flags )

! Note that expa only depends on temperature.
    tanh1 = tanh( ff / (( 2.0 * boltzmhz ) * t))
    DO i = 1, n_sps

      CALL create_beta ( catalog(i)%molecule, Catalog(i)%continuum, PB, T, &
        &  FF, Lines(Catalog(i)%Lines)%W, gl_slabs(n_ele,i), tanh1, bb )
      
      select case (catalog(i)%molecule)
      case (L_H2O)
        VMR = VMR_H2O
      case (L_O2)
        VMR = VMR_O2
      case (L_N2)
        VMR = VMR_N2
      case (L_O_18_O)
        VMR = VMR_O_18_O
      case (L_H2O_18)
        VMR = VMR_H2O_18
      case (L_O3)
        VMR = VMR_O3
      case (L_N2O)
        VMR = VMR_N2O
      case (L_HNO3)
        VMR = VMR_HNO3
      case default
        VMR=0.0_r8
      end select

      B = B + VMR*bb

    ENDDO

    ABSC=B/1000.0_r8     ! convert km-1 to m-1

    Call DestroyCompleteSlabs ( gl_slabs )
    Deallocate (pp, stat=status)
    if ( status /= 0 ) &
      & CALL MLSMessage(MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // ' pp ')
    Deallocate (tt, stat=status)
    if ( status /= 0 ) &
      & CALL MLSMessage(MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate // ' tt ')
    Deallocate (true_path_flags,stat=status)

  END SUBROUTINE get_beta_bill 
 
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

End Module Bill_GasAbsorption

! $Log$
! Revision 1.20  2003/07/04 02:49:12  vsnyder
! Simplify interface to Get_GL_Slabs_Arrays
!
! Revision 1.19  2003/06/18 14:43:26  bill
! modified interface to get_gl_slabs_arrays
!
! Revision 1.18  2003/05/16 23:53:36  livesey
! Now uses molecule indices rather than spectags
!
! Revision 1.17  2003/05/09 19:43:11  jonathan
! delete del_temp following Van's changes
!
! Revision 1.16  2003/05/06 20:45:45  jonathan
! some fixing after the merge
!
! Revision 1.15  2003/05/05 23:01:13  livesey
! Commented out call to create_beta for merge
!
! Revision 1.14  2003/02/11 00:48:22  jonathan
! change call to create_beta
!
! Revision 1.13  2003/02/06 00:21:10  jonathan
! change call to creat_beta
!
! Revision 1.12  2003/01/31 17:24:12  jonathan
! add Incl_Cld. cld_ext
!
! Revision 1.11  2003/01/30 18:25:06  pwagner
! Cosmetic changes; fixed bug where losVel was being changed
!
! Revision 1.10  2003/01/16 18:13:24  jonathan
! add Do_1D option to get_gl_slabs_arrays
!
! Revision 1.9  2002/12/13 02:07:18  vsnyder
! Use a SLABS structure for the slabs quantities, use spectag names
!
! Revision 1.8  2002/10/08 17:08:06  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.7  2002/08/22 00:13:01  jonathan
! upgrade to include more molecules
!
! Revision 1.6  2002/08/08 22:45:37  jonathan
! newly improved version
!
! Revision 1.5  2002/06/05 18:17:14  jonathan
!  fix bug
!
! Revision 1.4  2001/11/16 00:40:44  jonathan
! add losVel
!
! Revision 1.3  2001/11/15 00:55:02  jonathan
! fixed a bug
!
! Revision 1.2  2001/11/14 00:40:11  jonathan
! first working version
!
! Revision 1.1  2001/11/09 22:07:38  jonathan
! using Bill Read's new spec data
!




