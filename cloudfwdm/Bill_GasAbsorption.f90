! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Bill_GasAbsorption

  implicit NONE
  private
  public :: Get_Beta_Bill

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
  "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains

  subroutine Get_Beta_Bill (T, PB, F, RH, VMR_in, NS, ABSC, Catalog, LosVel )

!==============================================================
!      CALCULATE CLEAR-SKY ABSORPTION COEFFICIENT AT F AND T
!      USING BILL'S NEW SPECTROSCOPY DATA. 
!==============================================================

    use Get_Beta_Path_m, only: CREATE_BETA
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT, ALLOCATESLABS, &
      & DESTROYCOMPLETESLABS
    use MLSCommon, only: RK => R8, RP, IP
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
    REAL(rk), INTENT(IN) :: F               ! FREQUENCY IN GHz
    REAL(rk), INTENT(IN) :: T               ! TEMPERATURE (K)
    REAL(rk), INTENT(IN) :: PB              ! TOTAL AIR PRESSURE (hPa)
    REAL(rk), INTENT(IN) :: VMR_in(NS-1)    ! VMR OF SPECIES 
    REAL(rk), INTENT(IN) :: RH              ! H2O VOLUME MIXING RATIO
    REAL(rp), INTENT(IN) :: LosVel          !     OR RELATIVE HUMIDITY
    
    Type(Catalog_T), INTENT(IN) :: Catalog(:)
    Type (slabs_struct), DIMENSION(:,:), POINTER :: gl_slabs

    !------------------
    ! OUTPUTS
    !------------------
    REAL(rk), INTENT(out) :: ABSC           ! ABSORPTION COEFFICIENT (1/m)

    !------------------
    ! LOCAL VARIABLES
    !-----------------

    REAL(rk) :: B                           ! BETA (1/m/ppv)                 
    REAL(rk) :: FF                          ! FREQUENCY IN MHz               
    REAL(rp) :: myLosVel
    REAL(rk) :: P                           ! DRY AIR PARTIAL PRESSURE (hPa) 
    REAL(rk) :: VMR                         ! VOLUME MIXING RATIO            
    REAL(rk) :: VMR_H2O                     ! H2O VOLUME MIXING RATIO        
    REAL(rk) :: VMR_O3                      ! H2O VOLUME MIXING RATIO        
    REAL(rk) :: VMR_O2                      ! O2    VOLUME MIXING RATIO      
    REAL(rk) :: VMR_O_18_O                  ! O18O  VOLUME MIXING RATIO      
    REAL(rk) :: VMR_H2O_18                  ! H2O18 VOLUME MIXING RATIO      
    REAL(rk) :: VMR_N2                      ! N2 VOLUME MIXING RATIO         
    REAL(rk) :: VMR_N2O                     ! N2O VOLUME MIXING RATIO        
    REAL(rk) :: VMR_HNO3                    ! HNO3 VOLUME MIXING RATIO       
    REAL(rk) :: VP                          ! VAPOR PARTIAL PRESSURE (hPa)

    Integer(ip) :: n_sps, i
    real(rp) :: bb, tanh1
    real(rk), parameter :: Boltzmhz2 = 0.5 / h_over_k
    logical :: Do_1D

!-----------------------------------------------------------------------------
    if ( rh == 0.0_rk ) then
      p = pb
      vmr_h2o = vp / max(1.e-19_rk, p)
    else if ( rh == 110.0_rk ) then
      call RHtoEV ( T, rh, VP )
      P = PB - VP
      vmr_h2o = vp / max(1.e-19_rk, p)
    else
      vmr_h2o = rh                     ! RH HERE IS WATER VAPOR MIXING RATIO
      vp = vmr_h2o*pb                  ! VP IS VAPOR PRESSURE, PB IS TOTAL
      p = pb - vp                      ! PRESSURE, P IS DRY-AIR PRESSURE
    end if

    VMR_O2     = 0.2095_rk
    VMR_N2     = 0.805_rk
    VMR_O_18_O = VMR_O2*0.00409524_rk 
    VMR_H2O_18 = VMR_H2O*0.00204_rk
    VMR_O3     = VMR_in(1)
    VMR_N2O    = VMR_in(2)
    VMR_HNO3   = VMR_in(3)

    B = 0.0_rk
    FF = F*1000.0_rk

    n_sps = Size(Catalog)

    call allocateSlabs ( gl_slabs, 1, catalog, moduleName )

                              ! Bill uses km/sec 
    myLosVel=losVel*0.00_rp   ! The Doppler correction already been done 
                              ! in the FullCloudForwardModel, so set it 0

    call get_gl_slabs_arrays ( Catalog, (/ p /), (/ t /), myLosVel, gl_slabs, &
      & Do_1D, (/ .true. /) )

    ! Note that expa only depends on temperature.
    tanh1 = tanh( ff / ( boltzmhz2 * t))
    do i = 1, n_sps

      call create_beta ( catalog(i)%molecule, Catalog(i)%continuum, PB, T, &
        &  FF, Lines(Catalog(i)%Lines)%W, gl_slabs(1,i), tanh1, bb )
      
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
        VMR=0.0_rk
      end select

      B = B + VMR*bb

    end do

    ABSC = B/1000.0_rk     ! convert km-1 to m-1

    call DestroyCompleteSlabs ( gl_slabs )

  end subroutine get_beta_bill 

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Bill_GasAbsorption

! $Log$
! Revision 1.24  2003/07/09 23:30:16  vsnyder
! Interpret RH correctly, remove some old unused stuff
!
! Revision 1.23  2003/07/09 00:39:51  vsnyder
! Remove three unused variables
!
! Revision 1.22  2003/07/09 00:00:20  vsnyder
! Remove arrays known to have length 1
!
! Revision 1.21  2003/07/07 19:09:07  vsnyder
! Get Create_Beta from get_beta_path
!
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
