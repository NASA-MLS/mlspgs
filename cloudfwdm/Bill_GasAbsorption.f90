Module Bill_GasAbsorption

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

    use CREATE_BETA_M, only: CREATE_BETA
    use GLNP, only: NG
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT, ALLOCATEONESLABS, DESTROYCOMPLETESLABS
    use Make_Z_Grid_M, only: MAKE_Z_GRID
    use MLSCommon, only: R8, RP, IP
    use Molecules, only: SP_H2O, SP_H2O_18, SP_HNO3, SP_N2, SP_N2O, SP_O_18_O, &
      & SP_O2, SP_O3
    use SLABS_SW_M, only: GET_GL_SLABS_ARRAYS
    use SpectroscopyCatalog_m, only: CATALOG_T, LINES
    use WaterVapor, only: RHtoEV

    !-----------------
    ! INPUTS
    !-----------------

    INTEGER :: NS
    REAL(r8), INTENT(IN) :: F               ! FREQUENCY IN GHz
    REAL(r8) :: FF                          ! FREQUENCY IN MHz
    REAL(r8), INTENT(IN) :: T               ! TEMPERATURE (K)
    REAL(r8) :: P                           ! DRY AIR PARTIAL PRESSURE (hPa)
    REAL(r8) :: VP                          ! VAPOR PARTIAL PRESSURE (hPa)
    REAL(r8), INTENT(IN) :: PB              ! TOTAL AIR PRESSURE (hPa)
    REAL(r8) :: VMR_in(NS-1)                ! VMR OF SPECIES 
    REAL(r8), INTENT(IN) :: RH              ! H2O VOLUME MIXING RATIO OR RELATIVE HUMIDITY
    REAL(rp) :: LosVel
    
    Type(Catalog_T), INTENT(IN) :: Catalog(:)
    Type (slabs_struct), DIMENSION(:,:), POINTER :: gl_slabs

    !------------------
    ! OUTPUTS
    !------------------
    REAL(r8), INTENT(out) :: ABSC            ! ABSORPTION COEFFICIENT (1/m)

    !------------------
    ! LOCAL VARIABLES
    !-----------------

    REAL(r8) :: B                            ! BETA (1/m/ppv)
    REAL(r8) :: VMR                          ! VOLUME MIXING RATIO
    REAL(r8) :: VMR_H2O                      ! H2O VOLUME MIXING RATIO
    REAL(r8) :: VMR_O3                       ! H2O VOLUME MIXING RATIO
    REAL(r8) :: VMR_O2                       ! O2    VOLUME MIXING RATIO
    REAL(r8) :: VMR_O_18_O                   ! O18O  VOLUME MIXING RATIO
    REAL(r8) :: VMR_H2O_18                   ! H2O18 VOLUME MIXING RATIO
    REAL(r8) :: VMR_N2                       ! N2 VOLUME MIXING RATIO
    REAL(r8) :: VMR_N2O                      ! N2O VOLUME MIXING RATIO
    REAL(r8) :: VMR_HNO3                     ! HNO3 VOLUME MIXING RATIO

    Integer(ip) :: n_sps, n_path, i, j, k, m, nl, no_of_lines, n_ele, maxvert
    Integer(ip) :: Spectag, status
    REAL(rp) :: bb, v0, vm, tm, tp, bp, bm, del_temp
    REAL(rp), allocatable, dimension(:) :: PP, TT, z_psig
    logical :: Do_1D

!-----------------------------------------------------------------------------

    IF (RH .NE. 100.0_r8) THEN
       VMR_H2O = RH                     ! PH HERE IS WATER VAPOR MIXING RATIO
       VP=VMR_H2O*PB                    ! VP IS VAPOR PRESSURE, PB IS TOTAL
       P=PB-VP                          ! PRESSURE, P IS DRY-AIR PRESSURE
    ELSE IF(VMR_H2O .EQ. 100.0_r8) THEN
       CALL RHtoEV(T, 100.0_r8, VP)        ! RH HERE IS 100% RELATIVE HUMIDITY 
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

    B=0.0_r8
    FF = F*1000.0_r8

!    maxVert = Ng +1

    n_sps = Size(Catalog)
!    n_ele = 2*maxVert

    n_ele = 1  ! number of pressure levels along the path, in our case =1

    allocate ( gl_slabs (n_ele,n_sps), stat=status )
    allocate ( pp(n_ele), stat=status )
    allocate ( tt(n_ele), stat=status )

    do i = 1, n_sps
      no_of_lines =  size(Catalog(i)%Lines)
      gl_slabs(1:n_ele,i)%no_lines =  no_of_lines
      do j = 1, n_ele
          Call AllocateOneSlabs ( gl_slabs(j, i), no_of_lines )
      enddo
    enddo

    pp(1) = p
    tt(1) = t
    del_temp = 0.0_rp

    losVel=losVel*0.00_rp   ! Bill use km/sec ! The Doppler correction already been done 
                                              ! in the FullCloudForwardModel, so set it 0

    call get_gl_slabs_arrays( Catalog,PP(1:n_ele),TT(1:n_ele),losVel,gl_slabs,n_ele,del_temp, Do_1D)

    DO i = 1, n_sps
      Spectag = Catalog(i)%Spec_Tag

      CALL create_beta ( Spectag, Catalog(i)%continuum, PB, T, &
        &  FF, Lines(Catalog(i)%Lines)%W, gl_slabs(n_ele,i), bb )
      
      IF (Spectag .EQ. SP_H2O) THEN
        VMR = VMR_H2O
      ELSE IF (Spectag .EQ. SP_O2) THEN
        VMR = VMR_O2
      ELSE IF (Spectag .EQ. SP_N2) THEN
        VMR = VMR_N2
      ELSE IF (Spectag .EQ. SP_O_18_O) THEN
        VMR = VMR_O_18_O
      ELSE IF (Spectag .EQ. SP_H2O_18) THEN
        VMR = VMR_H2O_18
      ELSE IF (Spectag .EQ. SP_O3) THEN
        VMR = VMR_O3
      ELSE IF (Spectag .EQ. SP_N2O) THEN
        VMR = VMR_N2O
      ELSE IF (Spectag .EQ. SP_HNO3) THEN
        VMR = VMR_HNO3
      ELSE
        VMR=0.0_r8
      ENDIF

      B = B + VMR*bb

    ENDDO

    ABSC=B/1000.0_r8     ! convert km-1 to m-1

!    print*, B

    Call DestroyCompleteSlabs ( gl_slabs )
    DEAllocate(pp)
    DEAllocate(tt)

  END SUBROUTINE get_beta_bill 
 
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

End Module Bill_GasAbsorption

! $Log$
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




