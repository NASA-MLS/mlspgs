Module Bill_GasAbsorption

  use GLNP, only: NG
  use MLSCommon, only: r8, rp, ip
  use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT, ALLOCATEONESLABS
  use SpectroscopyCatalog_m, only: CATALOG_T, LINES
  use SLABS_SW_M, only: GET_GL_SLABS_ARRAYS
  use CREATE_BETA_M, only: CREATE_BETA
  use WaterVapor, only: RHtoEV
  Implicit NONE
  private
  PUBLIC :: get_beta_bill

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
  "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains

  SUBROUTINE get_beta_bill (T, PB, F, RH, VMR_O3, ABSC, Catalog, LosVel )

!==============================================================
!      CALCULATE CLEAR-SKY ABSORPTION COEFFICIENT AT F AND T
!      USING BILL'S NEW SPECTROSCOPY DATA. 
!      LATEST UPDATE: J.JIANG, NOVEMBER 14, 2001
!==============================================================

    !-----------------
    ! INPUTS
    !-----------------

    REAL(r8), INTENT(IN) :: F               ! FREQUENCY IN GHz
    REAL(r8) :: FF                          ! FREQUENCY IN MHz
    REAL(r8), INTENT(IN) :: T               ! TEMPERATURE (K)
    REAL(r8) :: P                           ! DRY AIR PARTIAL PRESSURE (hPa)
    REAL(r8) :: VP                          ! VAPOR PARTIAL PRESSURE (hPa)
    REAL(r8), INTENT(IN) :: PB              ! TOTAL AIR PRESSURE (hPa)
    REAL(r8) :: VMR_O3                      ! MINOR SPECIES 1-O3
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
    REAL(r8) :: VMR_O2                       ! O2    VOLUME MIXING RATIO
    REAL(r8) :: VMR_O_18_O                   ! O18O  VOLUME MIXING RATIO
    REAL(r8) :: VMR_H2O_18                   ! H2O18 VOLUME MIXING RATIO

    Integer(ip) :: n_sps, n_path, i, j, k, m, nl, no_of_lines, n_ele
    Integer(ip) :: Spectag, status
    REAL(rp) :: bb, v0, vm, tm, tp, bp, bm, del_temp
    REAL(rp), allocatable, dimension(:) :: LineWidth, PP, TT
    
!-----------------------------------------------------------------------------

    IF (RH .NE. 100._r8) THEN
       VMR_H2O = RH                     ! PH HERE IS WATER VAPOR MIXING RATIO
       VP=VMR_H2O*PB                    ! VP IS VAPOR PRESSURE, PB IS TOTAL
       P=PB-VP                          ! PRESSURE, P IS DRY-AIR PRESSURE
    ELSE IF(VMR_H2O .EQ. 100._r8) THEN
       CALL RHtoEV(T, 100._r8, VP)        ! RH HERE IS 100% RELATIVE HUMIDITY 
       P = PB-VP
       VMR_H2O = VP/(max(1.e-9_r8, P))
    END IF

!    VMR_H2O    = MAX(1.e-29_r8, VMR_H2O)
!    VMR_O3     = MAX(1.e-29_r8, VMR_O3)
    VMR_O2     = 0.209476_r8
    VMR_O_18_O = VMR_O2*0.00409524_r8 
    VMR_H2O_18 = VMR_H2O*0.00204_r8

    B=0._r8
    FF = F*1000._r8
    n_sps = Size(Catalog)
    n_ele = 1
    allocate ( gl_slabs (1,n_sps), stat=status )
    allocate ( pp(n_ele), stat=status )
    allocate ( tt(n_ele), stat=status )

    do i = 1, n_sps
      no_of_lines =  size(Catalog(i)%Lines)
      gl_slabs(1,i)%no_lines =  no_of_lines
      call AllocateOneSlabs ( gl_slabs(1, i), no_of_lines )
    enddo

    pp(1) = p
    tt(1) = t
    del_temp = 0.0_rp

    losVel=losVel*0.0001_rp   ! Bill use km/sec

    call get_gl_slabs_arrays(Catalog,PP(1:n_ele),TT(1:n_ele),losVel,gl_slabs,1,del_temp)


    DO i = 1, n_sps
      Spectag = Catalog(i)%Spec_Tag
      no_of_lines =  size(Catalog(i)%Lines)
      Allocate(LineWidth(no_of_lines))
      do k = 1, no_of_lines
        m = Catalog(i)%Lines(k)
        LineWidth(k) = Lines(m)%W
      end do

      CALL create_beta(Spectag, Catalog(i)%continuum, PB, T,                 &
        &  FF, no_of_lines, LineWidth, gl_slabs(1,i)%v0s, gl_slabs(1,i)%x1,  &
        &  gl_slabs(1,i)%y, gl_slabs(1,i)%yi, gl_slabs(1,i)%slabs1,  bb,     &
        &  gl_slabs(1,i)%dslabs1_dv0 )
      
      IF (Spectag .EQ. 18003) THEN
        VMR = VMR_H2O
      ELSE IF (Spectag .EQ. 32001) THEN
        VMR = VMR_O2
      ELSE IF (Spectag .EQ. 34001) THEN
        VMR = VMR_O_18_O
      ELSE IF (Spectag .EQ. 20003) THEN
        VMR = VMR_H2O_18
      ELSE IF (Spectag .EQ. 48004) THEN
        VMR = VMR_O3
      ENDIF

      B = B + VMR*bb

      DEAllocate(LineWidth)
    ENDDO

    ABSC=B/1000._r8     ! convert km-1 to m-1

!    print*, B

    DEAllocate(gl_slabs)
    DEAllocate(pp)
    DEAllocate(tt)

  END SUBROUTINE get_beta_bill 
 
End Module Bill_GasAbsorption

! $Log$
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




